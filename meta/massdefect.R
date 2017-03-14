#' Isotope extraction for single group of samples
#' 
#' 
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