install.packages('xlsx','rgl')
library(xlsx)
library(rgl)

# read in the data

data <- read.xlsx('./data/2D.xlsx',sheetIndex = 1)
raw <- data[,-1]
rownames(raw) <- data[,1]
rawm <- as.matrix(raw)

# get a matlab style
jet.colors <- 
    colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
# jet.colors <- colorRampPalette(c('yellow','red'))
# jet.colors <- colorRampPalette(c('blue','red'))
colorzjet <- jet.colors(1000)

# plot the surface map

rawm <- rawm
x <- 1:nrow(rawm)
z <- 100*(1:ncol(rawm))

open3d()
rgl.surface(x,z,rawm,color=colorzjet[findInterval(rawm, seq(min(rawm), max(rawm), length=1000))],back='lines')
axes3d('z',at=c(50,150,250,350,450,550,650,750,850,950,1050,1150),labels=c(4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75))
axes3d('x-+')
axes3d('y')
title3d('','','','response','pH')
mtext3d('time(second)','x-+',line=2)

# get the snapshot

rgl.snapshot("shot.png")
