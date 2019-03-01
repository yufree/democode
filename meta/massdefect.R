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
        # data <- as.data.frame(groups(xset))
        # z <- groupval(xset,'medret','into')
        # c <- apply(z,1,mean)
        # 
        df <- xcms::featureValues(xset, value = "into")
        c <- apply(df,1,mean)
        peaks <- xcms::featureDefinitions(xset)
        mz <- peaks$mzmed
        rt <- peaks$rtmed
        npeaks <- peaks$npeaks
        data <- cbind.data.frame(mzmed = mz, rtmed = rt, npeaks = npeaks)
        
        data$int <- c
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
#' Group peaks by the mass defect interval group for different substructures
#' @param list a peaks list with mass to charge
#' @param submass mass vector of sub structure of homologous series
#' @param mdgn mass defect groups numbers for interval, 20 means 0.05 inteval on mass defect scale from -0.5 to 0.5
#' @param lv group info for the data
#' @return list with mass defect analysis dataframe.
#' @export
getmdg <- function(list, submass = c(15.9949,14.003074,26.01568,14.01565,43.00581,30.01056,34.96885,78.91834), mdgn = 20, lv = NULL){
        if(is.null(list$stdmass)&is.null(list$paired)){
                mz <- list$mz
                rt <- list$rt
                data <- list$data
                colnames(data) <- lv
        }else if(is.null(list$stdmass)){
                mz <- list$mz[list$pairedindex]
                rt <- list$rt[list$pairedindex]
                data <- list$data[list$pairedindex,]
                colnames(data) <- lv
        }else{
                mz <- list$mz[list$stdmassindex]
                rt <- list$rt[list$stdmassindex]
                data <- list$data[list$stdmassindex,]
                colnames(data) <- lv
        }
        # perform mass defect analysis for std mass
        mda <- cbind.data.frame(mz = mz, rt = rt, data)
        for(i in 1:length(submass)){
                mdst <- round(submass[i])/submass[i]
                msdefect <- round(mz*mdst) - mz*mdst
                
                mdg <- cut(msdefect, seq(from = -.5, to = .5, by = 1/mdgn),include.lowest = T)
                mdg2 <- mdg[!is.na(mdg)]
                index <- mdg %in% enviGCMS::Mode(mdg2)
                
                name <- c(submass[i],paste0(submass[i],'gi'), 'majormdg')
                md <- cbind.data.frame(msdefect,mdg,index)
                colnames(md) <- name
                mda <- cbind.data.frame(mda,md)
        }
        # get the data
        list$mda <- mda
        return(list)
}
#' Paired correlationship among peak list based on cluster analysis
#' @param list a list with mzrt profile
#' @param rtcutoff cutoff of the distances in cluster
#' @param submass mass vector of sub structure of homologous series
#' @param mdcutoff mass defect cluster cutoff
#' @return list with retention time cluster, std mass defect analysis dataframe based on max average correlation
getcorstd <- function(list, rtcutoff = 9, submass = c(15.9949,14.003074,26.01568,14.01565,43.00581,30.01056,34.96885,78.91834), mdcutoff = 0.02){
        
        # paired mass diff analysis
        groups <- cbind.data.frame(mz = list$mz, rt = list$rt, list$data)
        resultstd <- resultsolo <- resultiso <- result <- NULL
        
        dis <- stats::dist(list$rt, method = "manhattan")
        fit <- stats::hclust(dis)
        rtcluster <- stats::cutree(fit, h=rtcutoff)
        n <- length(unique(rtcluster))
        message(paste(n, 'retention time cluster found.'))
        # search:
        for (i in 1:length(unique(rtcluster))) {
                # find the mass within RT
                rtxi <- list$rt[rtcluster == i]
                bin = groups[groups$rt %in% rtxi, ]
                medianrtxi <- stats::median(rtxi)
                
                if (nrow(bin) > 1) {
                        # get mz diff
                        cor <- stats::cor(t(bin[,-c(1,2)]))
                        cormean <- apply(cor,1,mean)
                        corindex <- which.max(cormean)
                        df <- cbind(bin[corindex,],rtg = i)
                        result <- rbind(result,df)
                }else{
                        solo <- cbind(bin,rtg = i)
                        resultsolo <- rbind(solo,resultsolo)
                }
        }
        
        resultstd <- rbind(result,resultsolo)
        resultstd <- unique(resultstd)
        
        # perform mass defect analysis for std mass
        for(i in 1:length(submass)){
                mdst <- round(submass[i])/submass[i]
                msdefect <- round(resultstd$mz*mdst) - resultstd$mz*mdst
                dis <- stats::dist(msdefect, method = "manhattan")
                fit <- stats::hclust(dis)
                mdcluster <- stats::cutree(fit, h=mdcutoff)
                n <- length(unique(mdcluster))
                message(paste(n, 'mass defect clusters found for mass', submass[i], 'substructures' ))
                name <- c(colnames(resultstd),submass[i],paste0(submass[i],'g'))
                resultstd <- cbind.data.frame(resultstd,msdefect,mdcluster)
                colnames(resultstd) <- name
        }
        
        # filter the list
        list$rtcluster <- rtcluster
        list$stdmassindex <- (round(list$mz,4) %in% round(resultstd$mz,4))
        list$stdmass <- resultstd
        return(list)
}
