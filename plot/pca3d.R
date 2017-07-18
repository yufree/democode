library(pca3d)
library(rgl)
library(readr)
library(car)

QDA_data_normalized <- read_csv("data/QDA_data_normalized.csv")
# QDA_data_normalized <- read_csv("data/QDA_data_processed.csv")

gr <- factor(t(QDA_data_normalized[,2]))
pca <- prcomp(QDA_data_normalized[,-c(1,2)],scale. = T,center = T)

groups <- gr

group.col <- c(rep("#FCB712",8),rep("#F37022",11),rep("#CD004D",9),rep("#6460AA",6),rep("#0089D1",5),rep("#0CB14B",8))
group.coln <- c("#FCB712","#F37022","#CD004D","#6460AA","#0089D1","#0CB14B")

# 1
pca3d(pca, col = group.col,group=groups,show.ellipses = T,axes.color = "black",show.group.labels = T,ellipse.ci = .95,show.plane = F, legend = 'topright')

# 2
pca3d <- function(pca,legend = F,groups, group.col = c(rep("#FCB712",8),rep("#F37022",11),rep("#CD004D",9),rep("#6460AA",6),rep("#0089D1",5),rep("#0CB14B",8)), group.coln = c("#FCB712","#F37022","#CD004D","#6460AA","#0089D1","#0CB14B")){
        levs <- levels(groups)
        x <- pca$x[,1]
        y <- pca$x[,2]
        z <- pca$x[,3]
        open3d()
        plot3d(x, y, z,xlab = 'PC 1 (23%)',ylab = 'PC 2 (18.3%)', zlab = 'PC 3 (13.7%)', col = group.col,type = 'p')
        for (i in 1:length(levs)) {
                group <- levs[i]
                selected <- groups == group
                xx <- x[selected]; yy <- y[selected]; zz <- z[selected]
                ellips <- ellipse3d(cov(cbind(xx,yy,zz)), 
                                    centre=c(mean(xx), mean(yy), mean(zz)), level = 0.95) 
                shade3d(ellips, col = group.coln[i], alpha = 0.1, lit = FALSE)
                # show group labels
                texts3d(mean(xx),mean(yy), mean(zz), text = group,col= group.coln[i], cex = 1)
        }
        if(legend){
                rgl::legend3d('topright',legend = levs,col = group.coln,pch = 19,bty = 'n',cex = 1)
        }
}
pca3d(pca,groups = gr)
#rgl.postscript("plot.pdf",fmt="pdf")

# 3
scatter3d(x=pca$x[,1],y=pca$x[,2],z = pca$x[,3],xlab = 'PC 1',ylab = 'PC 2', zlab = 'PC 3', col = group.col,ellipsoid = T,level = 0.95,surface = F,groups = gr)

## Final
library(pca3d)
library(rgl)
data_processed <- read_csv("C:/Users/SPMESUPERCPU/Downloads/data_processed.csv", col_types = cols(X1 = col_skip()))
col = c(rep("#FCB712",8),rep("#F37022",11),rep("#CD004D",9),rep("#6460AA",6),rep("#0089D1",5),rep("#0CB14B",8))
legend = 
        gr <- factor(t(data_processed[,1]))
pca <- prcomp(data_processed[,-1],scale. = T,center = T)
pca3d(pca, col = col,group=gr,show.ellipses = T,ellipse.ci = .75,show.plane = F,legend = 'topright')
snapshotPCA3d(file="fancy75.png")
rgl.postscript("graph.svg", fmt="svg")
pca3d(pca, col = col,group=gr,show.ellipses = T,axes.color = "black",show.group.labels = T,ellipse.ci = .95,show.plane = T
      , legend = 'topright')
snapshotPCA3d(file="fancy95.png")
pca3d(pca, group=gr,show.ellipses = T,show.group.labels = T,ellipse.ci = .50,show.plane = F, legend = 'topright')
snapshotPCA3d(file="fancy50.png")
rgl.postscript("graph.svg", fmt="svg")

