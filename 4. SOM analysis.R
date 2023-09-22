library(dplyr)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(sf)

setwd(".\\EVI")
########################## Normalise the data#############################

d.di.fill <- merge(DI.fill, eastmissi, by = c('CountryNum', 'StateNum', 'Route'))

d.di.fill <- d.di.fill %>% arrange(RouteID, Year) ###sort the data by routeid and year

TD <- d.di.fill[,c('RouteID', 'Year', 'TD')]
FD <- d.di.fill[,c('RouteID', 'Year', 'FD')]
PD <- d.di.fill[,c('RouteID', 'Year', 'PD')]

d.clusterid <- unique(TD$RouteID)

TD <- split(TD, TD$Year)
FD <- split(FD, FD$Year)
PD <- split(PD, PD$Year)


for (i in 1:length(TD)) {
  
  TD[[i]] <- TD[[i]][,3]
  
}

TD <- do.call(cbind, TD) 

colnames(TD) <- 1973:2016

for (i in 1:length(FD)) {
  
  FD[[i]] <- FD[[i]][,3]
  
}

FD <- do.call(cbind, FD) 

colnames(FD) <- 1973:2016

for (i in 1:length(PD)) {
  
  PD[[i]] <- PD[[i]][,3]
  
}

PD <- do.call(cbind, PD) 

colnames(PD) <- 1973:2016

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

TD.scale <- t(apply(TD, 1, range01))  
FD.scale <- t(apply(FD, 1, range01))
PD.scale <- t(apply(PD, 1, range01))

###########################################################################
library(kohonen)

##build a function to calculate topographic error

system('R CMD SHLIB ./src/map.c')###to load .c mapkohonen 
dyn.load("./map.dll")

topo.error.supersom <- function (somobj, type = c("nodedist", "bmu")) 
{
  if (attr(somobj,'class') != "kohonen") 
    stop("Topographic error measures as yet only defined for class 'kohonen'")
  type = match.arg(type)
  
  mdis <- sapply(1:length(somobj$data), function(i){
    
    dmat <- switch(type, nodedist = {
      as.matrix(dist(somobj$codes[[i]]))
    }, bmu = {
      if (!is.null(somobj$data)) data <- somobj$data[[i]]
      nd <- nrow(data)
      ncodes <- nrow(somobj$codes[[i]])
      np <- ncol(somobj$codes[[i]])
      
      distances <- .C("mapKohonen", as.double(data), as.double(somobj$codes[[i]]), 
                      as.integer(ncodes), as.integer(nd), as.integer(np), 
                      dists = double(nd * ncodes), NAOK = TRUE)$dists
      matrix(distances, nd, ncodes, byrow = TRUE)
    })
    
    dmat.ordered <- t(apply(dmat, 1, order))
    dists.2D <- unit.distances(somobj$grid)
    mean(mapply(function(i, j) dists.2D[i, j], dmat.ordered[,1], dmat.ordered[,2]))
  })
  
  return(mean(mdis))
  
}


set.seed(2019)

mygrid <- list(c(2,2),c(3,3),c(4,4),c(5,4),c(5,5))

############################calculate TE and QE#####################
tmpsom <- lapply(mygrid, function(grd) {
  
  set.seed(2019)
  
  TD.som <- som(list(TD.scale,FD.scale,PD.scale), grid = somgrid(grd[2], grd[1], topo = "hexagonal"), rlen = 1000,
                     normalizeDataLayers = F, dist.fcts = 'euclidean')
 
  ####togographic error 
  td.te <- c(topo.error.supersom(TD.som, 'bmu'), 'TD', 'TE', paste(grd[1], 'X', grd[2], sep = ' ' ))
  
  ###quantization error
  td.qe <- c(mean(TD.som$distances), 'TD', 'QE', paste(grd[1], 'X', grd[2], sep = ' ' ))
 
  dt <- rbind(td.te,td.qe)
  
  colnames(dt) <- c('value', 'diversity', 'err', 'samples')
  
  dt <- data.frame(dt)
  
  return(dt)
  
})

teqe <- data.table::rbindlist(tmpsom)

teqe$value <- as.numeric(levels(teqe$value))[teqe$value]

##plot the result of TE and QE

d.di <- teqe[teqe$diversity == 'TD',] 
  
d.te <- d.di[d.di$err == 'TE',]
  
d.qe <- d.di[d.di$err == 'QE',]
  
d.gg <- ggplot() +
    
    geom_line(size = 0.7, data = d.te, aes(y = value, x = samples, group = err, color = err)) +
    
    geom_point(size = 3, data = d.te, aes(y = value, x = samples, shape = err, color = err, fill = err)) +
    
    geom_text(aes(y = value, x = samples, label = scales::number(value, accuracy = 0.01)), data = d.te, vjust = -0.8) +
    
    geom_line(size = 0.7, data = d.qe, aes(y = value, x = samples, group = err, color = err)) +
    
    geom_point(size = 3, data = d.qe, aes(y = value, x = samples, shape = err, color = err, fill = err)) +
   
    geom_text(aes(y = value, x = samples, label = scales::number(value, accuracy = 0.01)), data = d.qe, vjust = 1.2) +
    
    geom_vline(xintercept = 2, linetype = 'dashed', colour = 'grey25') +
    
    scale_y_continuous(limits = c(1,1.6)) +
   
    theme_bw() +
    
    labs(x = 'Number of Neurons', y="Topographic/Quantization Error") +
    
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linetype = 3),
          axis.text.x = element_text(angle = 0, size = 8),
          axis.text.y = element_text(angle = 0, size = 10),
          legend.title = element_blank(),
          legend.position = c(0.85,0.9),
          legend.direction = "horizontal")



ggsave("qete.tiff", width = 9, height = 5, units = "in", plot = gg, dpi = 200)

########################################################################################
####################build a function to cluster and plot the som results################
########################################################################################

somcalplot <- function(grd){

set.seed(2019)

TD.som <- supersom(list(TD.scale,FD.scale,PD.scale), grid = somgrid(grd[2], grd[1], topo = "hexagonal"), rlen = 5000, 
                     normalizeDataLayers = F, dist.fcts = 'euclidean')
  
TDcluster_assignment <- TD.som$unit.classif

#plot the temporal trends of each pattern 

spdep(TDcluster_assignment, grd,BBS = T)

FPD <- d.di.fill[,c("RouteID","Year", "FD", "PD")]

FPD$FD <- do.call(c, lapply(unique(FPD$RouteID), function(x)range01(FPD[FPD$RouteID==x,]$FD)))
FPD$PD <- do.call(c, lapply(unique(FPD$RouteID), function(x)range01(FPD[FPD$RouteID==x,]$PD)))

TD <- data.frame(TD.scale)

TD$RouteID <- d.clusterid

##TD$RouteID <- retain

TD.list <- split(TD, TDcluster_assignment) ###group by the som results

##define function to calculate the stand error of mean

std1 <- function(x) mean(x) + sd(x)/sqrt(length(x))*1.96
std2 <- function(x) mean(x) - sd(x)/sqrt(length(x))*1.96

###define colours for the legend in ggplot

cols <- c("Taxonomic" = "#2e74b6", "Functional" = "#dd8039", "Phylogenetic" = "#f088bc", "Smoothed" = "grey50")

gg <- lapply(1:length(TD.list), function(x){
  
  tmp <- TD.list[[x]]
  
  tmpt <- tmp[,1:44] ##exclude the route id column
  
  tmpt <- split.default(tmpt, names(tmpt)) ##display by rows
  
  tmpt <- data.table::rbindlist(tmpt, idcol = T, use.names = F)
  
  names(tmpt) <- c("Year", "TD")
  
  tmpt$Year <- as.numeric(sub("^X", "", tmpt$Year)) ##unify the colnames
  
  tmpt$RouteID <- rep(tmp$RouteID, 44) ##add back the route id column
  
  tmpt <- merge(tmpt, FPD, by = c("Year", "RouteID"), all.x = T)
  
  tmp <- tmpt %>% group_by(Year) %>% 
    
    summarise(Rate = length(Year)/nrow(clustergeo), mTD = mean(TD), lTD = std1(TD), sTD = std2(TD),
              mFD = mean(FD), lFD = std1(FD), sFD = std2(FD), mPD = mean(PD), lPD = std1(PD), sPD = std2(PD))
  
  rate <- round(tmp$Rate[1],4)*100
  
  tmpt <- c(tmp$mFD, tmp$mPD, tmp$mTD)
  
  g <- factor(rep(c("FD", "PD", "TD"), each = nrow(tmp)))

  gg <- ggplot(tmp, aes(x = Year)) +
    
    stat_smooth(aes(y = mTD,  colour = 'Smoothed'), method = 'gam', formula = y ~ s(x)) +
    stat_smooth(aes(y = mFD,  colour = 'Smoothed'), method = 'gam', formula = y ~ s(x)) +
    stat_smooth(aes(y = mPD,  colour = 'Smoothed'), method = 'gam', formula = y ~ s(x)) +
    
    geom_line(aes(y = mTD, colour = "Taxonomic"), size = 0.8) +
    
    geom_errorbar(aes(ymin = sTD, ymax = lTD, colour = "Taxonomic"), width = 0.35) +
    
    geom_line(aes(y = mFD, colour = "Functional"), size = 0.8) +
    
    geom_errorbar(aes(ymin = sFD, ymax = lFD, colour = "Functional"), width = 0.35) +
    
    geom_line(aes(y = mPD, colour = "Phylogenetic"), size = 0.8) +
    
    geom_errorbar(aes(ymin = sPD, ymax = lPD, colour = "Phylogenetic"), width = 0.35) +
    
    scale_colour_manual(name="", values=cols, breaks = c("Taxonomic","Functional","Phylogenetic","Smoothed")) +
    
    #scale_fill_manual(name="", values=cols, breaks = c("Taxonomic","Functional","Phylogenetic","Smoothed")) +
    
    scale_x_continuous(expand = c(0.01,0), breaks = seq(1973, 2016, by = 10))+
    
    theme_bw() +
    
    labs(x = 'Year', y = "Diversity", subtitle = paste0("G",x,"; ",rate, "%: ", 100-rate,"%")) +
    scale_y_continuous(limits = c(0,1))+
    
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linetype = 3),
          axis.text.x = element_text(angle = 0, vjust = 0.6, size = 12),
          axis.text=element_text(colour = "black", size = 15),
          plot.subtitle = element_text(size = 14),
          legend.position = c(0.6,0.9),
          legend.title = element_text(size=16),
          legend.text = element_text(size = 16),
          legend.direction = "horizontal"
    )
  
  pgg <- gg + theme(legend.position = c(0.5,0.3))
  
  legendr <<- gtable_filter(ggplotGrob(pgg), "guide-box") 
  
  d.gg <- gg + theme(legend.position = 'none', axis.title = element_blank())
  
  return(d.gg)
  
})

Tdi <- grid.arrange(grobs = gg, ncol = grd[2], nrow = grd[1])

Tdi <- grid.arrange(Tdi, legendr,ncol = 1, heights = c(10,0.5))
                    #widths=unit.c(unit(1, "npc") - legendr$width, legendr$width))


ggsave(paste('nsem', grd[1], grd[2], '.tiff', sep = ''), width = grd[2]*4+2, height = grd[1]*4, units = "in", plot = Tdi, dpi = 200)# 6,8


}

################run the function
##
#########All BBS routes

usb <- read_sf("./NorAm.shp")

st_crs(usb) <- st_crs("ESRI:102003")

mygrid <- list(c(2,2),c(3,3),c(4,4),c(5,4),c(5,5))

for(p in 1:length(mygrid)){
  
  somcalplot(mygrid[[p]])
  
  cat(p,'\n')

  
}


####################routes located in the east of Mississipi

allpoints <- st_as_sf(clustergeo, coords = c("Longitude","Latitude"))

st_crs(allpoints) <- st_crs(4326)

mississipi <- sf::read_sf("./rivers.shp")

mississipi <- mississipi[mississipi$NAME=="Mississippi",1]
#st_crs(mississipi) <- st_crs(4326)

offsetX <- 50
offsetY <- 0

polyeast <- rbind(c(st_bbox(mississipi)['xmax'] + offsetX, 
                    st_bbox(mississipi)['ymin'] - offsetY),
                  c(st_bbox(mississipi)['xmax'] + offsetX, 
                    st_bbox(mississipi)['ymax'] + offsetY),
                  as.data.frame(st_coordinates(mississipi))[,c(1,2)],
                  c(st_bbox(mississipi)['xmax'] + offsetX, 
                    st_bbox(mississipi)['ymin'] - offsetY))  


polyeast <- st_sf("id" = 'sideEast', st_sfc(st_polygon(list(as.matrix(polyeast))), crs = 4326))

sf_use_s2(FALSE)
eastmissi <- st_intersection(allpoints, polyeast)

eastmissi$clusterID <- NULL

##################EVI

library(terra)

###########mask EVI to the North America
usb <- vect("./Noram.shp")

crs(usb) <- crs("ESRI:102003")

evi.list <- list.files(".\\Inputdata", pattern = ".tif$", full.names = T)

EVI <- rast(evi.list)

usb <- project(usb, crs(EVI))

usb <- crop(usb, ext(-179.9,0,0,90))

EVI<- mask(EVI,usb)

EVI <- crop(EVI, ext(usb))

names(EVI) <- 2000:2021

EVI <- EVI*0.0001

###############################preprare the dataframe
####all cells

lonlat <- xyFromCell(EVI, 1:(nrow(EVI)*ncol(EVI)))

evi.values <- extract(EVI,lonlat)

evi.values <- cbind(lonlat,evi.values)

evi.values$ID <- NULL
 
evi.values <- na.omit(evi.values)

lonlat <- data.frame(x=evi.values$x,y=evi.values$y)

evi.values$x <- evi.values$y <- NULL

evi.values <- as.matrix(evi.values)
#####cells overlap with BBS routes

allpoints <- vect(allpoints)

sub.evi.values <- extract(EVI,allpoints, ID=F, method = "bilinear")

sub.evi.values <- as.matrix(sub.evi.values)

coordspnts <- crds(allpoints,df=T)

#sub.evi.values <- t(apply(sub.evi.values, 1, range01))

######################SOM analysis

mygrid <- list(c(2,2),c(3,3),c(4,4),c(5,4),c(5,5))

usb <- read_sf("./NorAm.shp")

st_crs(usb) <- st_crs("ESRI:102003")

#####all cells

setwd("./AllCells")

for(p in 2:length(mygrid)){
  
  grd <- mygrid[[p]]
  
  set.seed(2019)
  
  evi.som <- supersom(evi.values, grid = somgrid(grd[2], grd[1], topo = "hexagonal"), rlen = 5000, 
                      normalizeDataLayers = F, dist.fcts = 'euclidean')
  
  evi.clusterID <- evi.som$unit.classif
  
  evi.df <- data.frame(Longitude = lonlat$x, Latitude = lonlat$y, clusterID = evi.clusterID)
  
  spdep(som.clusterID, grd, df = evi.df, BBS = F)
  
  cat(p,'\n')
  
  
}

###cells contain BBS routes

for(p in 1:length(mygrid)){
  
  grd <- mygrid[[p]]
  
  set.seed(2019)
  
  sub.evi.som <- supersom(sub.evi.values, grid = somgrid(grd[2], grd[1], topo = "hexagonal"), rlen = 5000, 
                      normalizeDataLayers = F, dist.fcts = 'euclidean')
  
  sub.evi.clusterID <- sub.evi.som$unit.classif
  
  sub.evi.df <- data.frame(Longitude = coordspnts$x, Latitude = coordspnts$y, clusterID = sub.evi.clusterID)
  
  spdep(som.clusterID, grd, df = sub.evi.df, BBS = F)
  
  cat(p,'\n')
  
  
}



 






















