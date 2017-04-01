getbrmd <- function(mass){
        sf <- 78/77.91051
        sm <- mass*sf
        sd <- ceiling(sm)-sm
        df <- as.data.frame(cbind(mass,sm,sd))
        plot(df$sd~df$sm,xlab = 'm/z',ylab = 'scaled MD')
        return(df)
}

getch2md <- function(mass){
        sf <- 14.00000/14.01565
        sm <- mass*sf
        sd <- ceiling(sm)-sm
        df <- as.data.frame(cbind(mass,sm,sd))
        plot(df$sd~df$sm,xlab = 'm/z',ylab = 'scaled MD')
        return(df)
}

getmassdefect <- function(mass,sf){
        sm <- mass*sf
        sd <- ceiling(sm)-sm
        df <- as.data.frame(cbind(mass,sm,sd))
        plot(df$sd~df$sm,xlab = 'm/z',ylab = 'scaled MD')
        return(df)
}

findbr <- function(xset,n=10,sf = 78/77.91051,step = 0.001,stepsd=0.003, cutoffint = 100, cutoffr=0.4){
        data <- as.data.frame(groups(xset))
        df <- peaks(xset)
        ind <- groupidx(xset)
        int <- NULL
        for(i in 1:length(ind)){
                int[i] <- mean(df[as.numeric(ind[[i]]),'into'])
        }
        data$int <- int
        data <- data[order(data$mzmed),]
        data <- data[data$npeaks >= n,]
        mz <- data$mzmed
        data$sm <- mz*sf
        data$sd <- ceiling(data$sm)-data$sm
        smstep <- seq(0,1,step)
        plot(data$sd~data$sm,type = 'n',xlab = 'm/z',ylab = 'scaled MD')
        result <- NULL
        for(i in 1:length(smstep)){
                mini = smstep[i]-stepsd
                maxi = smstep[i]+stepsd
                index = data$sd<maxi & data$sd>mini
                li <- data[index&data$int>cutoffint,]
                if(nrow(li)>2){
                        ratio <- max(li$int[li$int!=max(li$int)])/max(li$int)
                        if(ratio>cutoffr){
                                points(li$sd~li$sm,pch = 19,col = 'red')
                        }
                }else{
                        points(li$sd~li$sm,pch = 19,col = 'red')
                }
                result <- as.data.frame(rbind(result,li))
        }
        return(result[!duplicated(result), ])
}
#' Isotope extraction for single group of samples
getmassdiff <-
        function(xcmsSet,
                 massdiff,
                 rtwindow,
                 mzwindow,
                 ppm
        ) {
                # get group infomation
                groups = data.frame(xcmsSet@groups)
                peakIntensities = groupval(xcmsSet, method = "medret", value = intChoice)
                if (intChoice == "intb"){
                        peakIntensities[is.na(peakIntensities)] = 0
                }
                # order peaks by rt
                peakIntensities = peakIntensities[order(groups$rtmed), ]
                groups <- groups[order(groups$rtmed), ]
                groups$peakins <- apply(peakIntensities,1,mean)
                result <- NULL
                # search:
                for (i in 1:nrow(groups)) {
                        bin = groups[groups$rtmed-groups$rtmed[i] >= 0 & groups$rtmed-groups$rtmed[i] <= rtwindow, ]
                        if (nrow(bin) > 1) {
                                dis <- dist(bin$mzmed,method = 'manhattan')/massdiff
                                df <- data.frame(ms1=bin$mzmed[which(lower.tri(dis),arr.ind = T)[,1]],ms2 = bin$mzmed[which(lower.tri(dis),arr.ind = T)[,2]],diff=as.numeric(dis))
                                df$rdiff <- round(df$diff)
                                dfn <- df[df$diff<= df$rdiff * (1 + ppm / 1000000) + (df$ms1 * ppm /
                                                                                              1000000) / (massdiff * (1 - ppm / 1000000))
                                          &&
                                                  df$diff >= df$rdiff * (1 - ppm / 1000000) - (df$ms1 * ppm / 1000000) /
                                                  (massdiff * (1 + ppm / 1000000)),]
                                dfn$msdiff <- abs(dfn$ms1-dfn$ms2)
                                dfn <- dfn[dfn$msdiff < mzwindow,]
                                # candidate number of labeled atoms
                                result <- rbind(result,bin[bin$mzmed%in%dfn$ms1|bin$mzmed%in%dfn$ms2,])
                                result <- result[rownames(unique(result[,c("mzmed", "rtmed")])),]
                        }
                }
                result$sm <- result$mzmed*massdiff
                result$smd <- ceiling(result$sm)-result$sm
                return(result)
        }