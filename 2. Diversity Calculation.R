library(ape) 
library(picante)
library(FD)
library(ade4)
library(cluster)
library(vegan)
library(dplyr)

#####Prepare the trait data and build the functional dendrogram

###load the trait data and bird data then match them

Otraits <- openxlsx::read.xlsx('./BirdFun.xlsx')

traits <- Otraits[,c(8:19,24:31,36)]

traits$`Diet-Ver` <- rowSums(traits[,4:7])

traits[,4:7] <- NULL

traits <- traits[!is.na(traits$Scientific),]

######calculate the pairwise distances between using the Gower distance.

##fuzzy variable as the proportion

tabF <- traits[,c(3:8,18,9:15)]

tabFp <- prep.fuzzy(tabF, c(7,7), labels = c("Diet", "ForStrat"))

###quantity and dichotomous variable

tabQ <- data.frame(bodymass = traits[,17], row.names = rownames(traits)) ##body mass

tabD <- data.frame(pelagic = traits[,16], row.names = rownames(traits)) ##pelagic or not


###calculate the distance dissmissy

ktabl <- ktab.list.df(list(tabFp, tabQ, tabD))

distrait <- dist.ktab(ktabl, c("F", "Q", "D"))

######build the dendrogram using UPGMA clustering approach.

dendro <- hclust(distrait, method = "average")

dendro.ftree <- as.phylo(dendro)

######prepare the phylogenetic data

phylotree <- read.tree("./AllBirdsHackett1.tre")

######calculate the Diversity using Rao's Distance and sum of the branch lengths

###prepare the bird data, selecting data only including the routes satisfying gap <= 2 and 1969/73-2016

bird <- readRDS("./allinfo.rds")

bird <- merge(bird, specieslist, by = "AOU")

names(bird)[3:7] <- sub(pattern = ".x$", replacement = "", names(bird)[3:7])

bird[,c(15:17, 26:27)] <- NULL

route <- read.csv("./routeratain.csv")

route$X <- NULL

tmp <- anti_join(bird, route, by = c("CountryNum", "StateNum", "Route"))##return rows in bird not matched routes

bird <- setdiff(bird, tmp)##return routes of bird not in tmp, i.e. will be in routes

#### calculate the Diversity

routid <- unique(bird$RouteDataID)

pb <- txtProgressBar(1,length(routid), style = 3)

DI.list <- sapply(1:length(routid), function(r){
  
           setTxtProgressBar(pb,r)
  
           options(warn = -1)
  
           df <- bird[bird$RouteDataID == routid[r], ]
      
           CtryNum <- unique(df$CountryNum)
      
           StNum <- unique(df$StateNum)
      
           Rt <- unique(df$Route)
      
           yearlabel <- unique(df$Year)
      
           df <- merge(df,specieslist, by = "Spanish_Common_Name", all = T)
  
           df <- data.frame(df$SpeciesTotal, row.names = df$Spanish_Common_Name)
      
           df[,1][is.na(df[,1])] <- 0
      
           df <- t(df)
      
           df <- as.data.frame(df)
           
           phylo.c <- rep(0:100)
           
           phylo.uc <- rep(0,100)
           
           for (i in 1:100) {
             
             tmp <- sample(phylotree, 1)[[1]]
             
             tmp$`tip.label` <- sub("_"," ", tmp$`tip.label`)
             
             tmp <- keep.tip(tmp, tmp$`tip.label`[tmp$`tip.label` %in% specieslist$Spanish_Common_Name])
             
             phylo.c[i] <- raoD(df, tmp)$Dkk
             
             phylo.uc[i] <- pd(df,tmp)$PD
             
           }
           
           DI <- data.frame(TD = raoD(df)$Dkk, FD = raoD(df, dendro.ftree)$Dkk, PD = mean(phylo.c),
                            uTD = length(df[df!=0]), uFD = treedive(df, dendro, match.force = F), uPD = mean(phylo.uc),
                             CountryNum = CtryNum, StateNum = StNum, Route = Rt, Year = yearlabel)
      
           return(DI)
      
         }, simplify = F)  


DI <- bind_rows(DI.list)

write.csv(DI, "DIoriginal.csv")

DI <- read.csv("DIoriginal.csv")

### assign a unique ID for each route
 
DI <- DI %>% mutate(RouteID = group_indices_(DI, .dots = c("CountryNum", "StateNum", "Route")))

DI <- DI[with(DI, order(RouteID, Year)),]

### calculate the mean value for some route that were investigated more than one times

DI <- DI %>% group_by(RouteID, Year) %>% 
  summarise(TD = mean(TD), FD = mean(FD), PD = mean(PD), uTD = mean(uTD), uFD = mean(uFD), uPD = mean(uPD),
            CountryNum = unique(CountryNum), StateNum = unique(StateNum), Route = unique(Route))


###fill the data gap with NA

DI <- DI %>% 
  
  mutate(grp = cumsum(lag(Year > lead(Year, default = last(Year)), default = T))) %>%
  
  tidyr::complete(Year = min(Year):max(Year), grp) %>%
  
  arrange(grp) %>% 
  
  select(-grp)

DI <- DI[DI$Year >= 1973 & DI$Year <= 2016,]

write.csv(DI, 'diversity.csv')

#####distribution of the NAs

library(ggplot2)

nasum <- DI %>% group_by(Year) %>% summarise(na = sum(is.na(TD)))

ggplot(data = nasum, aes(x = Year, y = na)) +
  
  geom_bar(stat = 'identity') +
  
  scale_x_continuous(expand = c(0,0.2), breaks = seq(1973, 2016, by = 1)) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0, 15), breaks = seq(0, 15, by = 5))+
  
  theme_bw() +
  
  labs(x = 'Year', y = 'Nubmer of missing values') +
  
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = 3),
        axis.text.x = element_text(angle = 70, vjust = 0.6, size = 8),
        axis.text=element_text(colour = "black", size = 10)    )

ggsave("Nadistribution.png", width = 8, height = 4, units = "in", dpi = 200)

##### impute the missing data.

DI$TD <- as.numeric(DI$TD)
DI$FD <- as.numeric(DI$FD)
DI$PD <- as.numeric(DI$PD)
DI$uTD <- as.numeric(DI$uTD)
DI$uFD <- as.numeric(DI$uFD)
DI$uPD <- as.numeric(DI$uPD)

routid <- unique(DI$RouteID)

library(mice)

DIun <- subset(DI, select = -c(TD, FD, PD))

tmp.rf <- lapply(routid, function(x) {
  
  df <- DIun[DIun$RouteID == x,]
  
  df$CountryNum <- df$CountryNum[!is.na(df$CountryNum)][1]
  df$StateNum <- df$StateNum[!is.na(df$StateNum)][1]
  df$Route <- df$Route[!is.na(df$Route)][1]
  
  tmpt <- mice(df, m = 50, method = "rf", maxit = 50, printFlag = F)
  
  tmpt <- complete(tmpt)
  
  x <<- x
  
  cat(x,'\n')
  
  return(tmpt)
  
})

DI.fill <- bind_rows(tmp.rf)  

write.csv(DI.fill,"./DIfill_norm.csv")


