data <- read.csv('data/GY.csv')
dm <- as.matrix(data[,c('PFBA','PFPeA',"PFHxA","PFOA","PFBS","PFHxS",'PFOS')])
heatmap(dm)
dmlog <- log(dm+1)
rownames(dmlog) <- data[,1]
heatmap(dmlog)
par(mar=c(4,4,2,2))
library('pheatmap')
pheatmap(dmlog,color = colorRampPalette(c("white",'navy'))(100),fontsize = 18)
