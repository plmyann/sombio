
############## BBS data preprocessing
##1. summarise the annual number of routes
##2. the spatial distribution of the total number of routes in each adminstration region
##3. retain routes including year from 19?? to 2017


library(dplyr)
library(data.table)
library(ggplot2)
library(rgdal)
library(gganimate)

input <- './DataInput'

setwd(input)

###read all BBS data

allcsv <- list.files(input,'.csv$')

tmp <- sapply(allcsv, read.csv, simplify = F)

alldf <- rbindlist(tmp)

##########summarise count of observations by year and route ID

alldf.sum <- alldf %>% 
  
  group_by(RouteDataID, Year) %>%
  
  count() %>% 
  
  group_by(Year) %>% 
  
  count()

######merge the information weather and route to each observation

routes <- read.csv('./routes.csv')

weather <- read.csv('./weather.csv')

rouwther <- merge(routes, weather, by = c('CountryNum', 'StateNum', 'Route'), all = T)

alldfallinfo <- merge(alldf, rouwther, by = 'RouteDataID', all = T)

alldfallinfo <- alldfallinfo[!is.na(RouteDataID)]

saveRDS(alldfallinfo, 'allinfo.rds')

routyear <- alldfallinfo %>% group_by(CountryNum.x, StateNum.x, Year.x) %>% summarise(n = n_distinct(RouteDataID))

colnames(routyear) <- c("CountryNum", "StateNum", "Year", "n")

region <- read.csv('./RegionCode.csv')

routyear <- merge(routyear,region, by = c("CountryNum", "StateNum"))

write.csv(routyear,'routyear.csv')

#########################plot the spatial distribution of routes

dates <- paste(routyear$Year, "01", "01", sep = "-")

dates_unique <- unique(dates)[order(unique(dates))]

##usb1 <- us_boundaries(dates[1])

usb1 <- readOGR("./NorAm.shp")

usb1 <- spTransform(usb1, CRS("ESRI:102003"))#Albers 1983

sel <- usb1@data$ADMIN_NAME %in% routyear$StateName

summary(sel)

usb1$ADMIN_NAME[!sel]

usbj <- left_join(usb1@data, routyear, by = c('ADMIN_NAME' = 'StateName'))

for(i in 1:length(dates_unique)) {
  
  usbi <- usb1
  
  print(st_crs(usbi))
  
  routyearj <- routyear[routyear$Year == lubridate::year(dates_unique[i]),]
  
  usbi@data = left_join(usb1@data, routyearj, by = c('ADMIN_NAME' = 'StateName'))
  
  usbi.map <- fortify(usbi, region = 'ADMIN_NAME')
  
  usbi.map <- merge(usbi.map, usbi@data, by.x = 'id', by.y = 'ADMIN_NAME')
  
  usbj <- plyr::rbind.fill(usbj, usbi.map)
}


usbj <- usbj[usbj$long <= 0,]

usb1966 <- usbj[usbj$Year <= 1970,]

 p <- ggplot(data = usbj, aes(x = long, y = lat, group = group, fill = n)) + 
  
        geom_polygon(color = 'grey') +
  
        labs(title = 'Year: {current_frame}', x = 'Longitude', y = 'Latitude', fill = 'Number of routes') +
  
        theme_bw() + 
   
        scale_fill_gradient(high = "#132B43", low = "#56B1F7") +
  
        transition_manual(Year)
 
 animate(p, fps = 1, renderer = gifski_renderer('anrout1.gif', width = 800*5, height = 1500*5))
 
 ### build a funciton to select consecutive numbers
 
 consecutive <- function(t){
   
   t <- sort(as.numeric(unique(t)))
   
   b <- cumsum(c(1, diff(t)>1))
  
   c <- rle(b)
   
   return(c)
 }
 
####################################
 ##### retain data including 1969-2016 with gap <= 2.
 
 route.sum <- list()
 
 for (y in 1966:1975) {
   
   weathery <- weather[weather$Year >= y,]
   
   tmp <- weathery %>%
     
     group_by(CountryNum, StateNum, Route) %>%
     
     summarise(Nyear = n_distinct(Year), Minyear = min(Year), Maxyear = max(Year), 
               MinyearID = RouteDataID[which.min(Year)], MaxyearID = RouteDataID[which.max(Year)],
               gap = length(consecutive(Year)$lengths)-1)
   ###gap is the number of interruptions, Ngap is the number of gap years
   
   route.sum[[y-1965]] <- data.frame(tmp)
   
 }

 route.sum <- rbindlist(route.sum) 
 
 route.sum$Ngap <-  route.sum$Maxyear - route.sum$Minyear - route.sum$Nyear + 1
 
 ####plot distribution in different starting year and gaps 
 
 tmp <- route.sum[route.sum$Minyear <= 1975 & route.sum$Maxyear >= 2016,]
 
 df <- tmp %>% group_by(Minyear, Ngap) %>% count()
 
 df <- df %>% group_by(Minyear) %>% mutate(nn = cumsum(n))
 
 ymax <- max(df$n)+50
 
 df$Minyear <- as.factor(df$Minyear)
 
 gg1 <- ggplot(df, aes(Ngap, n, fill = Minyear, color = Minyear)) +
   geom_line(size = 0.6)+
   scale_x_continuous(expand = c(0,0), breaks = seq(0, max(df$Ngap), by = 1))+
   scale_y_continuous(expand = c(0,0), limits = c(0, ymax), 
                      breaks = seq(0, ymax, by = 50))+
   theme_bw() + 
   labs(x = 'Number of gaps', y = 'Nubmer of routes', fill = 'Start Year', color = 'Start Year') +
   theme(panel.border = element_blank(),
         axis.line = element_line(colour = "black",size =0.5),
         plot.margin=unit(c(0.5,0.5,0.2,0.2), "cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         axis.title.x = element_text(colour = "black",size = 8, margin=margin(10,0,0,0)),
         axis.title.y = element_text(colour = "black",size = 12, margin=margin(0,10,0,0)),
         #axis.title.y = element_blank(),
         axis.text=element_text(colour = "black", size = 7),
         axis.ticks=element_line(size = 0.5),
         plot.title = element_text(size=12,margin=margin(0,10,10,10), hjust = 0.5)) 
 
 ggsave("gapsum.tiff", width = 8, height = 4, units = "in", plot = gg1, dpi = 200)
 
 
##### retain data including 1969-2016 with gap <= 2.

routeretain1 <- tmp[tmp$Minyear <= 1973 & tmp$Minyear >= 1969 & tmp$Maxyear >= 2016 & tmp$Ngap <= 2,]

routeretain <- unique(routeretain1, by = c("CountryNum","StateNum","Route"))

write.csv(routeretain, './routeratain.csv')
 
