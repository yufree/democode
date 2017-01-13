# devtools::install_github('yufree/sva-devel')
library(xcms)
library(RColorBrewer)
library(sva)
library(limma)
library(CAMERA)
library(qvalue)

svacor <- function(xset,lv,annotation=F,polarity = "positive",nSlaves=12){
        data <- groupval(xset,"maxint", value='into')
        mod <- model.matrix(~lv)
        mod0 <- as.matrix(c(rep(1,ncol(data))))
        svafit <- sva(data,mod)
        if (svafit$n.sv == 0){
            svaX <- model.matrix(~lv)
            lmfit <- lmFit(data,svaX)
            Signal<-lmfit$coef[,1:nlevels(lv)]%*%t(svaX[,1:nlevels(lv)])
            error <- data-Signal
            rownames(Signal) <- rownames(error) <- rownames(data)
            colnames(Signal) <- colnames(error) <- colnames(data)
            pValues = f.pvalue(data,mod,mod0)
            qValues = qvalue(pValuesSv)
            qValues = qValues$qvalues
            if(annotation){
                dreport <- annotateDiffreport(xset,metlin = T,polarity = polarity,nSlaves = nSlaves)
                dreport <- dreport[order(as.numeric(rownames(dreport))),]
                return(list(data,Signal,error,pValues,qValues,dreport))
            }else{
                return(list(data,Signal,error,pValues,qValues))
            }
        }else{
            message('data is correcting ...')
            svaX <- model.matrix(~lv+svafit$sv)
            lmfit <- lmFit(data,svaX)
            Batch<- lmfit$coef[,(nlevels(lv)+1):(nlevels(lv)+svafit$n.sv)]%*%t(svaX[,(nlevels(lv)+1):(nlevels(lv)+svafit$n.sv)])
            Signal<-lmfit$coef[,1:nlevels(lv)]%*%t(svaX[,1:nlevels(lv)])
            error <- data-Signal-Batch
            datacor <- Signal+error
            rownames(datacor) <- rownames(Signal) <- rownames(Batch) <- rownames(error) <- rownames(data)
            colnames(datacor) <- colnames(Signal) <- colnames(Batch) <- colnames(error) <- colnames(data)
            modSv = cbind(mod,svafit$sv)
            mod0Sv = cbind(mod0,svafit$sv)
            pValuesSv = f.pvalue(data,modSv,mod0Sv)
            qValuesSv = qvalue(pValuesSv)
            qValuesSv = qValuesSv$qvalues
            
            pValues = f.pvalue(data,mod,mod0)
            qValues = qvalue(pValues)
            qValues = qValues$qvalues
            if(annotation){
                dreport <- annotateDiffreport(xset,metlin = T,polarity = polarity,nSlaves = nSlaves)
                dreport <- dreport[order(as.numeric(rownames(dreport))),]
                return(list(data,datacor,Signal,Batch,error,pValues,qValues,pValuesSv,qValuesSv,dreport,svafit$pprob.gam,svafit$pprob.b))
            }else{
                return(list(data,datacor,Signal,Batch,error,pValues,qValues,pValuesSv,qValuesSv,svafit$pprob.gam,svafit$pprob.b))
            }
        }     
}

svapca <- function(list){
        data <- list[[1]]
        Signal <- list[[2]]
        Batch <- list[[3]]
        error <- list[[4]]
        
        par(mfrow=c(2,4),mar = c(2.75, 2.2, 2.6, 1))
        pcao <- prcomp(t(data), center=TRUE, scale=TRUE)
        plot(pcao, type = "l",main = "PCA")
        
        pca <- prcomp(t(Signal), center=TRUE, scale=TRUE) 
        plot(pca, type = "l",main = "PCA-signal")
        
        pcab <- prcomp(t(Batch), center=TRUE, scale=TRUE)
        plot(pcab, type = "l",main = "PCA-batch")
        
        pcae <- prcomp(t(error), center=TRUE, scale=TRUE)
        plot(pcae, type = "l",main = "PCA-error")
        
        plot(pcao$x[,1], 
             pcao$x[,2], 
             xlab="PC1",
             ylab="PC2",
             pch=colnames(data),
             cex=2,
             main = "PCA")
        
        plot(pca$x[,1], 
             pca$x[,2], 
             xlab="PC1",
             ylab="PC2",
             pch=colnames(Signal),
             cex=2,
             main = "PCA-signal")
        
        plot(pcab$x[,1], 
             pcab$x[,2], 
             xlab="PC1",
             ylab="PC2",
             pch=colnames(Batch),
             cex=2,
             main = "PCA-batch")
        
        plot(pcae$x[,1], 
             pcae$x[,2], 
             xlab="PC1",
             ylab="PC2",
             pch=colnames(error),
             cex=2,
             main = "PCA-error")
}

svaplot <- function(list, pqvalues="sv"){
        data <- list[[1]]
        Signal <- list[[2]]
        Batch <- list[[3]]
        error <- list[[4]]
        pValues <- list[[5]]
        qValues <- list[[6]]
        pValuesSv <- list[[7]]
        qValuesSv <- list[[8]]
        par(mfrow=c(1,5),mar = c(3,3,1,1))
        icolors <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)
        
        if(pqvalues == "anova" & sum(pValues<0.05&qValues<0.05)!=0){
                zlim <- range(c(Signal[pValues<0.05&qValues<0.05,],data[pValues<0.05&qValues<0.05,],Batch[pValues<0.05&qValues<0.05,],error[pValues<0.05&qValues<0.05,]))
    
                image(t(data[pValues<0.05&qValues<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/(ncol(data)-1)), labels=colnames(data),cex.axis=0.618,las=2)
    
                image(t(Signal[pValues<0.05&qValues<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks-signal',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/(ncol(data)-1)), labels=colnames(Signal),cex.axis=0.618,las=2)
    
                image(t(Batch[pValues<0.05&qValues<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks-batch',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/(ncol(data)-1)), labels=colnames(Batch),cex.axis=0.618,las=2)
    
                image(t(error[pValues<0.05&qValues<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks-error',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/(ncol(data)-1)), labels=colnames(error),cex.axis=0.618,las=2)
    
                breaks <- seq(zlim[1], zlim[2], round((zlim[2]-zlim[1])/10))
                poly <- vector(mode="list", length(icolors))
                plot(1,1,t="n",xlim=c(0,1), ylim=zlim,xaxt='n', yaxt='n',xaxs="i", yaxs="i",ylab = '',xlab = 'intensity',frame.plot=F)
                axis(4,at=breaks,labels = round(breaks),las=1,pos=0.4)
                bks <- seq(zlim[1], zlim[2], length.out=(length(icolors)+1))
                for(i in seq(poly)){
                        polygon(c(0.1,0.1,0.3,0.3), c(bks[i], bks[i+1], bks[i+1], bks[i]),  col=icolors[i], border=NA)
                        }
                return(list(Signal[pValues<0.05&qValues<0.05,],pValues<0.05&qValues<0.05))
                }
        else if(pqvalues == "sv"&sum(pValuesSv<0.05&qValuesSv<0.05)!=0){
                zlim <- range(c(Signal[pValuesSv<0.05&qValuesSv<0.05,],data[pValuesSv<0.05&qValuesSv<0.05,],Batch[pValuesSv<0.05&qValuesSv<0.05,],error[pValuesSv<0.05&qValuesSv<0.05,]))

                image(t(data[pValuesSv<0.05&qValuesSv<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/(ncol(data)-1)), labels=colnames(data),cex.axis=0.618,las=2)

                image(t(Signal[pValuesSv<0.05&qValuesSv<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks-signal',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/(ncol(data)-1)), labels=colnames(Signal),cex.axis=0.618,las=2)

                image(t(Batch[pValuesSv<0.05&qValuesSv<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks-batch',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/(ncol(data)-1)), labels=colnames(Batch),cex.axis=0.618,las=2)

                image(t(error[pValuesSv<0.05&qValuesSv<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks-error',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/(ncol(data)-1)), labels=colnames(error),cex.axis=0.618,las=2)

                breaks <- seq(zlim[1], zlim[2], round((zlim[2]-zlim[1])/10))
                poly <- vector(mode="list", length(icolors))
                plot(1,1,t="n",xlim=c(0,1), ylim=zlim,xaxt='n', yaxt='n',xaxs="i", yaxs="i",ylab = '',xlab = 'intensity',frame.plot=F)
                axis(4,at=breaks,labels = round(breaks),las=1,pos=0.4)
                bks <- seq(zlim[1], zlim[2], length.out=(length(icolors)+1))
                for(i in seq(poly)){
                        polygon(c(0.1,0.1,0.3,0.3), c(bks[i], bks[i+1], bks[i+1], bks[i]),  col=icolors[i], border=NA)
                }
                return(list(Signal[pValuesSv<0.05&qValuesSv<0.05,],pValuesSv<0.05&qValuesSv<0.05))
                }
        else{
                zlim <- range(c(Signal,data,Batch,error))

                image(t(data),col=icolors,xlab = 'samples',ylab = 'peaks',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/(ncol(data)-1)), labels=colnames(data),cex.axis=0.618,las=2)
    
                image(t(Signal),col=icolors,xlab = 'samples',ylab = 'peaks-signal',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/(ncol(data)-1)), labels=colnames(Signal),cex.axis=0.618,las=2)
    
                image(t(Batch),col=icolors,xlab = 'samples',ylab = 'peaks-batch',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/(ncol(data)-1)), labels=colnames(Batch),cex.axis=0.618,las=2)
    
                image(t(error),col=icolors,xlab = 'samples',ylab = 'peaks-error',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/(ncol(data)-1)), labels=colnames(error),cex.axis=0.618,las=2)
    
                breaks <- seq(zlim[1], zlim[2], round((zlim[2]-zlim[1])/10))
                poly <- vector(mode="list", length(icolors))
                plot(1,1,t="n",xlim=c(0,1), ylim=zlim,xaxt='n', yaxt='n',xaxs="i", yaxs="i",ylab = '',xlab = 'intensity',frame.plot=F)
                axis(4,at=breaks,labels = round(breaks),las=1,pos=0.4)
                bks <- seq(zlim[1], zlim[2], length.out=(length(icolors)+1))
                for(i in seq(poly)){
                        polygon(c(0.1,0.1,0.3,0.3), c(bks[i], bks[i+1], bks[i+1], bks[i]),  col=icolors[i], border=NA)
                }
        }
}
