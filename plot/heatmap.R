data <- read.csv('data/GY.csv')
dm <- as.matrix(data[,c('PFBA','PFPeA',"PFHxA","PFOA","PFBS","PFHxS",'PFOS')])
heatmap(dm)
dmlog <- log(dm+1)
rownames(dmlog) <- data[,1]
heatmap(dmlog)
par(mar=c(4,4,2,2))
library('pheatmap')
pheatmap(dmlog,color = colorRampPalette(c("white",'navy'))(100),fontsize = 18)

x <- c(0,0,0,0,0,0,0,0,0,0,0,-1,-2,-3,1,2,3,3,3,-3,-3,-2,-2,2,2)
y <- c(0,1,2,3,4,5,-1,-2,-3,-4,-5,0,0,0,0,0,0,-5,5,5,-5,-3,3,3,-3)
z <- c(0,0,2600000,520000,280000,2730,930000,450000,410000,0,0,0,0,0,0,0,1100000,0,0,0,0,320,350000,440000,480000)
z1 <- c(0,0,310000,1100000,290000,13000,310000,900000,530000,0,0,0,0,0,0,0,740000,0,0,0,0,1710,430000,530000,520000)
z2 <- c(0,0,1700000,1400000,550000,5080,880000,1200000,960000,0,0,0,0,0,0,0,1600000,0,0,0,0,1470,730000,850000,840000)
z3 <- c(0,0,1900000,3300000,1300000,32000,1200000,3200000,2100000,0,0,0,0,0,0,0,3400000,0,0,0,0,19000,1600000,1800000,2100000)
z4 <- c(0,0,190000,68000,36000,1120,63000,60000,51000,0,0,0,0,0,0,0,150000,0,0,0,0,230,37000,48000,45000)
z5 <- c(0,0,32000,43000,19000,690,18000,47000,31000,0,0,0,0,0,0,0,55000,0,0,0,0,350,24000,28000,30000)
z6 <- c(0,0,32000,43000,19000,690,18000,47000,31000,0,0,0,0,0,0,0,55000,0,0,0,0,350,24000,28000,30000)

library('fields')
png('Salbutamol.png',width = 600,height = 500)
quilt.plot(x,y,z,nx=15,ny=15,cex = 3,main = 'Salbutamol')
dev.off()

png('Sertaline.png',width = 600,height = 500)
quilt.plot(x,y,z1,nx=15,ny=15,cex = 3,main = 'Sertaline')
dev.off()

png('Cocaine.png',width = 600,height = 500)
quilt.plot(x,y,z2,nx=15,ny=15,cex = 3,main = 'Cocaine')
dev.off()

png('Fentanyl.png',width = 600,height = 500)
quilt.plot(x,y,z3,nx=15,ny=15,cex = 3,main = 'Fentanyl')
dev.off()

png('Oxycodone.png',width = 600,height = 500)
quilt.plot(x,y,z4,nx=15,ny=15,cex = 3,main = 'Oxycodone')
dev.off()

png('Buprenorphine.png',width = 600,height = 500)
quilt.plot(x,y,z5,nx=15,ny=15,cex = 3,main = 'Buprenorphine')
dev.off()

png('Clenbuterol.png',width = 600,height = 500)
quilt.plot(x,y,z6,nx=15,ny=15,cex = 3,main = 'Clenbuterol')
dev.off()

