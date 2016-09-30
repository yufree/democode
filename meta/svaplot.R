library(multtest)
library(xcms)
library(CAMERA)
# library(faahKO)
# library(rafalib)
library(RColorBrewer)
library(sva)
library(limma)

svadata <- function(path,...){
        cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
        xset <- xcmsSet(cdffiles,nSlaves = 12)
        xset <- group(xset)
        xset2 <- retcor(xset, method = "obiwarp")
        # you need group the peaks again for this corrected data
        xset2 <- group(xset2)
        xset3 <- fillPeaks(xset2,nSlaves=12)
        data <- groupval(xset3,"maxint", value='into',...)
        return(data)
}

svaplot <- function(data,lv,pqvalues=F){
        mod <- model.matrix(~lv)
        mod0 <- as.matrix(c(rep(1,ncol(data))))
        svafit <- sva(data,mod)
        svaX <- model.matrix(~lv+svafit$sv)
        lmfit <- lmFit(data,svaX)
        
        Batch<- lmfit$coef[,3:(2+svafit$n.sv)]%*%t(svaX[,3:(2+svafit$n.sv)])
        Signal<-lmfit$coef[,1:2]%*%t(svaX[,1:2])
        error <- data-Signal-Batch
        rownames(Signal) <- rownames(Batch) <- rownames(error) <- rownames(data)
        colnames(Signal) <- colnames(Batch) <- colnames(error) <- colnames(data)
        modSv = cbind(mod,svafit$sv)
        mod0Sv = cbind(mod0,svafit$sv)
        pValuesSv = f.pvalue(data,modSv,mod0Sv)
        qValuesSv = p.adjust(pValuesSv,method="BH")
        
        pValues = f.pvalue(data,mod,mod0)
        qValues = p.adjust(pValues,method = "BH")
        
        dataout <- cbind(data,pValues,qValues,pValuesSv,qValuesSv)
        par(mfrow=c(1,5),mar = c(2.75, 2.2, 2.6, 1))
        icolors <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)
        if(pqvalues == "anova"){
                zlim <- range(c(Signal[pValues<0.05&qValues<0.05,],data[pValues<0.05&qValues<0.05,],Batch[pValues<0.05&qValues<0.05,],error[pValues<0.05&qValues<0.05,]))
                
                image(t(data[pValues<0.05&qValues<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/11), labels=colnames(data),cex.axis=0.8,las=2)
                
                image(t(Signal[pValues<0.05&qValues<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks-signal',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/11), labels=colnames(Signal),cex.axis=0.8,las=2)
                
                image(t(Batch[pValues<0.05&qValues<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks-batch',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/11), labels=colnames(Batch),cex.axis=0.8,las=2)
                
                image(t(error[pValues<0.05&qValues<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks-error',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/11), labels=colnames(error),cex.axis=0.8,las=2)
                
                breaks <- seq(zlim[1], zlim[2], round((zlim[2]-zlim[1])/10))
                poly <- vector(mode="list", length(icolors))
                plot(1,1,t="n",xlim=c(0,1), ylim=zlim,xaxt='n', yaxt='n',xaxs="i", yaxs="i",ylab = '',xlab = 'intensity',frame.plot=F)
                axis(4,at=breaks,labels = round(breaks),las=1,pos=0.4)
                bks <- seq(zlim[1], zlim[2], length.out=(length(icolors)+1))
                for(i in seq(poly)){
                        polygon(c(0.1,0.1,0.3,0.3), c(bks[i], bks[i+1], bks[i+1], bks[i]),  col=icolors[i], border=NA)
                }
                return(dataout[pValues<0.05&qValues<0.05,])
        }else if(pqvalues == "sv"){
                zlim <- range(c(Signal[pValuesSv<0.05&qValuesSv<0.05,],data[pValuesSv<0.05&qValuesSv<0.05,],Batch[pValuesSv<0.05&qValuesSv<0.05,],error[pValuesSv<0.05&qValuesSv<0.05,]))
                
                image(t(data[pValuesSv<0.05&qValuesSv<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/11), labels=colnames(data),cex.axis=0.8,las=2)
                
                image(t(Signal[pValuesSv<0.05&qValuesSv<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks-signal',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/11), labels=colnames(Signal),cex.axis=0.8,las=2)
                
                image(t(Batch[pValuesSv<0.05&qValuesSv<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks-batch',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/11), labels=colnames(Batch),cex.axis=0.8,las=2)
                
                image(t(error[pValuesSv<0.05&qValuesSv<0.05,]),col=icolors,xlab = 'samples',ylab = 'peaks-error',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/11), labels=colnames(error),cex.axis=0.8,las=2)
                
                breaks <- seq(zlim[1], zlim[2], round((zlim[2]-zlim[1])/10))
                poly <- vector(mode="list", length(icolors))
                plot(1,1,t="n",xlim=c(0,1), ylim=zlim,xaxt='n', yaxt='n',xaxs="i", yaxs="i",ylab = '',xlab = 'intensity',frame.plot=F)
                axis(4,at=breaks,labels = round(breaks),las=1,pos=0.4)
                bks <- seq(zlim[1], zlim[2], length.out=(length(icolors)+1))
                for(i in seq(poly)){
                        polygon(c(0.1,0.1,0.3,0.3), c(bks[i], bks[i+1], bks[i+1], bks[i]),  col=icolors[i], border=NA)
                }
                return(dataout[pValuesSv<0.05&qValuesSv<0.05,])
        }else{
                zlim <- range(c(Signal,data,Batch,error))
                
                image(t(data),col=icolors,xlab = 'samples',ylab = 'peaks',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/11), labels=colnames(data),cex.axis=0.8,las=2)
                
                image(t(Signal),col=icolors,xlab = 'samples',ylab = 'peaks-signal',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/11), labels=colnames(Signal),cex.axis=0.8,las=2)
                
                image(t(Batch),col=icolors,xlab = 'samples',ylab = 'peaks-batch',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/11), labels=colnames(Batch),cex.axis=0.8,las=2)
                
                image(t(error),col=icolors,xlab = 'samples',ylab = 'peaks-error',xaxt="n",yaxt="n",zlim=zlim)
                axis(1, at=seq(0,1,1/11), labels=colnames(error),cex.axis=0.8,las=2)
                
                breaks <- seq(zlim[1], zlim[2], round((zlim[2]-zlim[1])/10))
                poly <- vector(mode="list", length(icolors))
                plot(1,1,t="n",xlim=c(0,1), ylim=zlim,xaxt='n', yaxt='n',xaxs="i", yaxs="i",ylab = '',xlab = 'intensity',frame.plot=F)
                axis(4,at=breaks,labels = round(breaks),las=1,pos=0.4)
                bks <- seq(zlim[1], zlim[2], length.out=(length(icolors)+1))
                for(i in seq(poly)){
                        polygon(c(0.1,0.1,0.3,0.3), c(bks[i], bks[i+1], bks[i+1], bks[i]),  col=icolors[i], border=NA)
                }
                return(dataout)
        }
}

