getRAW <- function(xset,lv){
        lsa <- list()
        index <- 0
        m <- length(table(lv))
        for(i in 1:m){
                if(i == 1){
                        index <- index
                }else{
                        index <- index+table(lv)[i-1]
                }
                ls <- list()
                for(j in 1:table(lv)[i]){
                        RAW <- getXcmsRaw(xset,sampleidx=j+index)
                        RAWdata <- RAW@env$profile
                        colnames(RAWdata) <- RAW@scantime
                        rownames(RAWdata) <- seq(RAW@mzrange[1],RAW@mzrange[2],0.1)
                        ls[[j]] <- RAWdata
                }
                lsa[[i]] <- ls
        }
        return(lsa)
}

getalign <- function(list){
        a <- lapply(list,rownames)
        c <- Reduce(intersect,a)
        e <- lapply(list,function(x){x[c,]})
        return(e)
}

getQCraw <- function(path,mzrange,rtrange,index=NULL){
	cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
	if(index){
          cdffiles <- cdffiles[index]
      }
    nsamples <- length(cdffiles)
    area <- numeric()
    for (i in 1:nsamples){
    	RAW <- xcmsRaw(cdffiles[i])
    	peak <- rawEIC(RAW,mzrange,rtrange)
    	area[i] <- sum(peak$intensity)
    }
    return(area)
}

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
