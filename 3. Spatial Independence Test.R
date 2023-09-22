library(rgdal)
library(maptools)
library(spatstat)
library(dplyr)
library(ggplot2)

##build a function to test spatial independence
spdep <- function(cluster_assignment, grd, df, BBS=T){
  
northamr.proj <- CRS("ESRI:102003")##

if(BBS==T){
  
  routdf <- read.csv('./routes.csv')##all routes
  
  dindex <- data.frame(cluster_assignment)
  
  clusterid <- data.frame(clusterID = dindex[,1], RouteID = d.clusterid)
  
  clusterdf <- left_join(clusterid, d.di.fill, by = 'RouteID')
  
  clusterdf <- clusterdf %>% group_by(clusterID,CountryNum, StateNum, Route) %>% count()
  
  clusterdf$n <- NULL
  
  clustergeo <- merge(clusterdf, routdf, by=c('CountryNum','StateNum', 'Route'), all.y = F)
  
  
} else{
  
  clustergeo <- data.frame(df)
}


clustergeo$clusterID <- as.factor(clustergeo$clusterID)

d.p <- rep(NA, n_distinct(clustergeo$clusterID))

gg <- lapply(1:n_distinct(clustergeo$clusterID), function(clu){
  
  d.geo <- clustergeo[clustergeo$clusterID == clu,]
  f.geo <- clustergeo[clustergeo$clusterID != clu,]
  
##############plot the density histogram
  d.spnts <- SpatialPoints(cbind(d.geo$Longitude,d.geo$Latitude), proj4string = CRS("+proj=longlat +datum=WGS84"))
  f.spnts <- SpatialPoints(cbind(f.geo$Longitude,f.geo$Latitude), proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  d.spnts.proj <- spTransform(d.spnts, northamr.proj)
  f.spnts.proj <- spTransform(f.spnts, northamr.proj)
  
  d.pnts <- as(d.spnts.proj, 'ppp')
  
  d.pnts.sf <- st_as_sf(d.spnts.proj)
  f.pnts.sf <- st_as_sf(f.spnts.proj)
  
  if(d.pnts$n < 5){
    
    d.p[clu] <- NA
    
    ##############plot the routes distribution
    
    qq <- ggplot() +
      
      geom_sf(data = usb, color = 'grey80', fill = "#9ddac9") +
      
      labs(x = 'Longitude', y = 'Latitude', title = paste0("G",clu,"; ",round(nrow(d.geo)/nrow(clustergeo),4)*100,'% : ', 
                                                           (1-round(nrow(d.geo)/nrow(clustergeo),4))*100, "%",
                                                           '; P = ', d.p[clu])) +
      
      #geom_sf(shape = 21, color = "white",fill = "#dd8039", size = 3, data = f.pnts.sf, alpha = 0.5) +
      
      #geom_sf(shape = 21, color = "white",fill = "#2e74b6", size = 3, data = d.pnts.sf, alpha = 0.8)+
      
      geom_sf(shape = 16, color =  "#dd8039", size = .1, data = f.pnts.sf, alpha = 0.5) +
      
      geom_sf(shape = 16, color = "#2e74b6", size = .1, data = d.pnts.sf, alpha = 0.8)+
      
      
      scale_x_continuous(breaks = c(-70,-90,-110,-130)) + 
      
      theme(plot.title = element_text(size=16),
            axis.text = element_text(size = 12))+
      
      theme_bw() 
    
    ggsave(paste0('sdd',grd[1], grd[2], '_', clu, '.tiff'), width = 10, height = 10, units = "in", plot = qq, dpi = 200)
    
    d.qq <- qq + theme(legend.position = 'none', axis.title = element_blank())
    
    
  }else if(d.pnts$n >= 5) {
  
  d.poly <- as(adehabitatHR::mcp(d.spnts.proj, unout = 'm', percent = 100), 'owin')
  
  marks(d.pnts) <- NULL
  
  Window(d.pnts) <- d.poly
  
  d.ann <- mean(nndist(d.pnts, k=1))
  
  n <- 1000
  
  d.ann.r <- vector(length = n)
  
  set.seed(2019)
  
  for (i in 1:n) {
    
    d.ran.p <- rpoint(n = d.pnts$n, win = d.poly)
    
    d.ann.r[i] <- mean(nndist(d.ran.p, k = 1))
    
  }
  
  #plot(d.ran.p, pch = 16, main = NULL, cols = rgb(0,0,0,.5))
  
  n.greater <- sum(d.ann.r > d.ann)
  
  p <- min(n.greater + 1, n+1-n.greater)/(n+1)*2
  
  d.p[clu] <- scales::number(p, accuracy = 0.001)
  
  d.ann.r <- data.frame(val = d.ann.r/1000)
  
  rangednn <- range(d.ann/1000, d.ann.r$val)
  
  d.gg <- ggplot(data = d.ann.r, aes(x = val)) +
    
    geom_density(alpha = 0.5, adjust = 1, fill = "#f088bc",color = "#f088bc") +
     
    scale_x_continuous(expand = c(0,0.02), limits = range(d.ann/1000, d.ann.r$val), 
                       breaks = c(round(min(rangednn),-1),round(sum(rangednn)/2,-1),round(max(rangednn),-1))) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(density(d.ann.r$val)$y)+0.005))+
     
    labs(x = 'ANN distances (km)', y = 'Frequency', subtitle =paste('P: ', scales::number(p, accuracy = 0.001))) +
     
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank()) +
     
    geom_vline(xintercept = d.ann/1000, colour = 'gray20', size =0.8) +
    geom_vline(xintercept = mean(d.ann.r$val), colour = 'gray20', size =0.8, linetype="dashed") 
  
  ggsave(paste('sdep',grd[1], grd[2], '_', clu, '.tiff', sep = ''), width = 6, height = 6, units = "in", plot = d.gg, dpi = 200)
  
  d.gg <- d.gg + labs(x = 'ANN distances(km)', y = 'Frequency', subtitle = '') + 
    theme(axis.title=element_text(size=14),
          axis.text = element_text(size=10))
  
##############plot the routes distribution
  
  qq <- ggplot() +
    
    geom_sf(data = usb, color = 'grey80', fill = "#9ddac9") +
    
    labs(x = 'Longitude', y = 'Latitude', title = paste0("G",clu,"; ",round(nrow(d.geo)/nrow(clustergeo),4)*100,'% : ', 
                                                         (1-round(nrow(d.geo)/nrow(clustergeo),4))*100, "%",
                                                         '; P = ', d.p[clu])) +
    
    #geom_sf(shape = 21, color = "white",fill = "#dd8039", size = 3, data = f.pnts.sf, alpha = 0.5) +
    
    #geom_sf(shape = 21, color = "white",fill = "#2e74b6", size = 3, data = d.pnts.sf, alpha = 0.8)+
    
    geom_sf(shape = 16, color =  "#dd8039", size = .1, data = f.pnts.sf, alpha = 0.5) +
    
    geom_sf(shape = 16, color = "#2e74b6", size = .1, data = d.pnts.sf, alpha = 0.8)+
    
    scale_x_continuous(breaks = c(-70,-90,-110,-130)) + 
    
    theme_bw() +
    
    theme(plot.title = element_text(size=16),
          axis.text = element_text(size = 12))
  
  qqq <- ggplot()+annotation_custom(ggplotGrob(qq), xmin = -.05, xmax = 1.05 , ymin = -.065, ymax = 1.07) + 
    annotation_custom(ggplotGrob(d.gg), xmin =.1, xmax = .4 , ymin = .06, ymax = .5) + theme_minimal()+
    theme(plot.margin=unit(c(0,0,0,0), "cm"))
    #theme(plot.margin =  unit(c(-1.6, -1.2, -1.5, -1.4), "cm"))
  
  ggsave(paste('sdd',grd[1], grd[2], '_', clu, '.tiff', sep = ''), width = 10, height = 8, units = "in", plot = qqq, dpi = 200)
  
  #d.qq <- qqq + theme(legend.position = 'none', axis.title = element_blank())
  d.qq <- ggplot()+annotation_custom(ggplotGrob(qq + theme(axis.title = element_blank())), xmin = -.05, xmax = 1.05 , ymin = -.057, ymax = 1.07) + 
    annotation_custom(ggplotGrob(d.gg), xmin =.1, xmax = .4 , ymin = .06, ymax = .5) + theme_minimal()+
    theme(plot.margin=unit(c(0,0,0,0), "cm"))
  
  }
  
  return(d.qq) 
  
})

Tdi <- grid.arrange(grobs =gg, ncol = grd[2], nrow = grd[1])

Tdi <- grid.arrange(Tdi, ncol = 1, heights = unit.c(unit(0.96, "npc")),
                    #widths=unit.c(unit(1, "npc") - legendr$width, legendr$width),
                    bottom = textGrob("Longitude", gp = gpar(cex = 1.8), vjust = -0.2),
                    left = textGrob("Latitude", rot = 90, vjust = 0.8, gp = gpar(cex = 1.8)))

#Tdi <- cowplot::plot_grid(plotlist = lapply(gg, "+", theme(plot.margin=unit(c(0,0,0,0), "cm"))), ncol = grd[2], nrow = grd[1])

ggsave(paste0('allspdep', grd[1], grd[2], '.tiff'), width = grd[2]*5, height = grd[1]*4, units = "in", plot = Tdi, dpi = 300)# 6,8



}

