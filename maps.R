data <- read.table('clipboard',header=T,sep="\t")
install.packages(c('maptools','raster'))
library(maptools)
library(raster)
adm <- getData('GADM', country='China', level=3)
mar<-(adm[adm$NAME_3=="Shouguang",])
adm <- getData('SRTM', lon=118, lat=37)
plot(adm)
points(data$E,data$N,pch=19)
plot(mar, axes=T,xlim=c(118.6,119),ylim=c(37.1,37.3),add=T)

# China maps

adm <- getData('GADM', country='China', level=1)
adm1 <- getData('GADM', country='Taiwan', level=1)
plot(adm)
plot(adm1, add = TRUE)
