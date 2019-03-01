#' Get xcmsset object in one step with optimized methods.
#' @param path the path to your data
#' @param index the index of the files
#' @param BPPARAM used for BiocParallel package
#' @param pmethod parameters used for different instrumentals such as 'hplcorbitrap', 'uplcorbitrap', 'hplcqtof', 'hplchqtof', 'uplcqtof', 'uplchqtof'. The parameters were from the references
#' @param minfrac minimum fraction of samples necessary in at least one of the sample groups for it to be a valid group, default 0.67
#' @param ... arguments for xcmsSet function
#' @details the parameters are extracted from the papers. If you use name other than the name above, you will use the default setting of XCMS. Also I suggest IPO packages or apLCMS packages to get reasonable data for your own instrumental. If you want to summit the results to a paper, remember to include those parameters.
#' @return a xcmsset object for that path or selected samples
#' @references Patti, G. J.; Tautenhahn, R.; Siuzdak, G. Nat. Protocols 2012, 7 (3), 508–516.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' }
#' @export
getdata <-
        function(path,
                 index = F,
                 BPPARAM = BiocParallel::SnowParam(),
                 pmethod = "hplcorbitrap",
                 minfrac = 0.67,
                 ...) {
                cdffiles <- list.files(path, recursive = TRUE,
                                       full.names = TRUE)
                if (index) {
                        cdffiles <- cdffiles[index]
                }
                if (pmethod == "hplcorbitrap") {
                        xset <- xcms::xcmsSet(
                                cdffiles,
                                BPPARAM = BPPARAM,
                                method = "centWave",
                                ppm = 2.5,
                                peakwidth = c(10,
                                              60),
                                prefilter = c(3, 5000),
                                ...
                        )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else {
                                xset <- xcms::group(
                                        xset,
                                        bw = 5,
                                        mzwid = 0.015,
                                        minfrac = min
                                )
                                xset2 <- xcms::retcor(xset, 'obiwarp')
                                # you need group the peaks again for this corrected
                                # data
                                xset2 <-
                                        xcms::group(
                                                xset2,
                                                bw = 5,
                                                mzwid = 0.015,
                                                minfrac = minfrac
                                        )
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == "uplcorbitrap") {
                        xset <- xcms::xcmsSet(
                                cdffiles,
                                BPPARAM = BPPARAM,
                                method = "centWave",
                                ppm = 2.5,
                                peakwidth = c(5,
                                              20),
                                prefilter = c(3, 5000),
                                ...
                        )
                        xset <-
                                xcms::group(
                                        xset,
                                        bw = 2,
                                        mzwid = 0.015,
                                        minfrac = minfrac
                                )
                        xset2 <- xcms::retcor(xset, 'obiwarp')
                        # you need group the peaks again for this corrected
                        # data
                        xset2 <-
                                xcms::group(
                                        xset2,
                                        bw = 2,
                                        mzwid = 0.015,
                                        minfrac = minfrac
                                )
                        xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                } else if (pmethod == "hplcqtof") {
                        xset <- xcms::xcmsSet(
                                cdffiles,
                                BPPARAM = BPPARAM,
                                method = "centWave",
                                ppm = 30,
                                peakwidth = c(10,
                                              60),
                                prefilter = c(0, 0),
                                ...
                        )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else {
                                xset <- xcms::group(
                                        xset,
                                        bw = 5,
                                        mzwid = 0.025,
                                        minfrac = minfrac
                                )
                                xset2 <- xcms::retcor(xset, 'obiwarp')
                                # you need group the peaks again for this corrected
                                # data
                                xset2 <-
                                        xcms::group(
                                                xset2,
                                                bw = 5,
                                                mzwid = 0.025,
                                                minfrac = minfrac
                                        )
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == "hplchqtof") {
                        xset <- xcms::xcmsSet(
                                cdffiles,
                                BPPARAM = BPPARAM,
                                method = "centWave",
                                ppm = 15,
                                peakwidth = c(10,
                                              60),
                                prefilter = c(0, 0),
                                ...
                        )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else {
                                xset <- xcms::group(
                                        xset,
                                        bw = 5,
                                        mzwid = 0.015,
                                        minfrac = minfrac
                                )
                                xset2 <- xcms::retcor(xset, 'obiwarp')
                                # you need group the peaks again for this corrected
                                # data
                                xset2 <-
                                        xcms::group(
                                                xset2,
                                                bw = 5,
                                                mzwid = 0.015,
                                                minfrac = minfrac
                                        )
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == "uplcqtof") {
                        xset <- xcms::xcmsSet(
                                cdffiles,
                                BPPARAM = BPPARAM,
                                method = "centWave",
                                ppm = 30,
                                peakwidth = c(5,
                                              20),
                                prefilter = c(0, 0),
                                ...
                        )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else {
                                xset <- xcms::group(
                                        xset,
                                        bw = 2,
                                        mzwid = 0.025,
                                        minfrac = minfrac
                                )
                                xset2 <- xcms::retcor(xset, 'obiwarp')
                                # you need group the peaks again for this corrected
                                # data
                                xset2 <-
                                        xcms::group(
                                                xset2,
                                                bw = 2,
                                                mzwid = 0.025,
                                                minfrac = minfrac
                                        )
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == "uplchqtof") {
                        xset <- xcms::xcmsSet(
                                cdffiles,
                                BPPARAM = BPPARAM,
                                method = "centWave",
                                ppm = 15,
                                peakwidth = c(5,
                                              20),
                                prefilter = c(0, 0),
                                ...
                        )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else {
                                xset <- xcms::group(
                                        xset,
                                        bw = 2,
                                        mzwid = 0.015,
                                        minfrac = minfrac
                                )
                                xset2 <- xcms::retcor(xset, 'obiwarp')
                                # you need group the peaks again for this corrected
                                # data
                                xset2 <-
                                        xcms::group(
                                                xset2,
                                                bw = 2,
                                                mzwid = 0.015,
                                                minfrac = minfrac
                                        )
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else {
                        xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM,
                                              ...)
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else {
                                xset <- xcms::group(xset, minfrac = minfrac)
                                xset2 <- xcms::retcor(xset, 'obiwarp')
                                # you need group the peaks again for this corrected
                                # data
                                xset2 <-
                                        xcms::group(xset2, minfrac = minfrac)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                }
                return(xset3)
        }

#' Get the csv files to be submitted to Metaboanalyst
#' @param xset the xcmsset object which you want to submitted to Metaboanalyst
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param name file name
#' @return dataframe with data needed for Metaboanalyst if your want to perform local analysis.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' xset <- getdata(cdfpath, pmethod = ' ')
#' getupload(xset)
#' }
#' @export
getupload <- function(xset,
                      method = "medret",
                      intensity = "into",
                      name = "Peaklist") {
        peakIntensities <- xcms::groupval(xset, method,
                                          intensity)
        peakIntensities[is.na(peakIntensities)] = 0

        data <-
                rbind(group = as.character(xcms::phenoData(xset)$class),
                      peakIntensities)
        data <- data[!duplicated(rownames(data)),]
        filename <- paste0(name, ".csv")
        utils::write.csv(data, file = filename)
        return(data)
}

#' output the similarity of two dataset
#' @param xset1 the first dataset
#' @param xset2 the second dateset
#' @return similarity on retention time and rsd %
#' @export
getsim <- function(xset1, xset2) {
        data1 <- gettechrep(xset1)[, c("mzmed", "rtmed",
                                       "rsd")]
        data2 <- gettechrep(xset2)[, c("mzmed", "rtmed",
                                       "rsd")]

        # data1$weight <- ifelse(data1$rsd > 100, 0, 0.2)
        # data1$weight[data1$rsd < 80] <- 0.4
        # data1$weight[data1$rsd < 60] <- 0.6
        # data1$weight[data1$rsd < 40] <- 0.8
        # data1$weight[data1$rsd < 20] <- 1 data1$rtorder
        # <- order(data1$rtmed)
        data1$mzmedn <- round(data1$mzmed, 0.1)

        # data2$weight <- ifelse(data2$rsd > 100, 0, 0.2)
        # data2$weight[data2$rsd < 80] <- 0.4
        # data2$weight[data2$rsd < 60] <- 0.6
        # data2$weight[data2$rsd < 40] <- 0.8
        # data2$weight[data2$rsd < 20] <- 1 data2$rtorder
        # <- order(data2$rtmed)
        data2$mzmedn <- round(data2$mzmed, 0.1)

        data <- merge(data1, data2, by = "mzmedn")
        data <- data[stats::complete.cases(data),]
        cor1 <- cor(data$rtmed.x, data$rtmed.y)
        cor2 <- cor(data$rsd.x, data$rsd.y)
        cor <- c(cor1, cor2)
        return(cor)
}

#' get the data of QC compound for a group of data
#' @param path data path for your QC samples
#' @param mzrange mass of the QC compound
#' @param rtrange retention time of the QC compound
#' @param index index of the files contained QC compounds, default is all of the compounds
#' @return number vector, each number indicate the peak area of that mass and retention time range
#' @export
getQCraw <- function(path, mzrange, rtrange, index = NULL) {
        cdffiles <- list.files(path, recursive = TRUE,
                               full.names = TRUE)
        if (index) {
                cdffiles <- cdffiles[index]
        }
        nsamples <- length(cdffiles)
        area <- numeric()
        for (i in 1:nsamples) {
                RAW <- xcms::xcmsRaw(cdffiles[i])
                peak <- xcms::rawEIC(RAW, mzrange, rtrange)
                area[i] <- sum(peak$intensity)
        }
        return(area)
}

#' Isotope extraction for single group of samples with certain mass diff
#' @param xcmsSet  a xcmsSet object
#' @param massdiff mass defection
#' @param rtwindow retention time range
#' @param mzwindow mass charge ratio window
#' @param ppm resolution of the mass spectrum
#' @return table with mass, retention time, scaled mass and scaled mass defect
getmassdiff <- function(xcmsSet, massdiff, rtwindow,
                        mzwindow, ppm) {
        # get group infomation
        groups = data.frame(xcmsSet@groups)
        peakIntensities = xcms::groupval(xcmsSet, "medret",
                                         "inio")
        # order peaks by rt
        peakIntensities = peakIntensities[order(groups$rtmed),]
        groups <- groups[order(groups$rtmed),]
        groups$peakins <- apply(peakIntensities, 1, mean)
        result <- NULL
        # search:
        for (i in 1:nrow(groups)) {
                bin = groups[groups$rtmed - groups$rtmed[i] >=
                                     0 &
                                     groups$rtmed - groups$rtmed[i] <= rtwindow,]
                if (nrow(bin) > 1) {
                        dis <- stats::dist(bin$mzmed, method = "manhattan") / massdiff
                        df <-
                                data.frame(
                                        ms1 = bin$mzmed[which(lower.tri(dis),
                                                              arr.ind = T)[, 1]],
                                        ms2 = bin$mzmed[which(lower.tri(dis),
                                                              arr.ind = T)[, 2]],
                                        diff = as.numeric(dis)
                                )
                        df$rdiff <- round(df$diff)
                        dfn <-
                                df[df$diff <= df$rdiff * (1 + ppm / 1e+06) +
                                           (df$ms1 * ppm / 1e+06) / (massdiff * (1 -
                                                                                         ppm /
                                                                                         1e+06)) && df$diff >= df$rdiff *
                                           (1 - ppm / 1e+06) - (df$ms1 * ppm /
                                                                        1e+06) / (massdiff *
                                                                                          (1 + ppm /
                                                                                                   1e+06)),]
                        dfn$msdiff <- abs(dfn$ms1 - dfn$ms2)
                        dfn <- dfn[dfn$msdiff < mzwindow,]
                        # candidate number of labeled atoms
                        result <- rbind(result, bin[bin$mzmed %in%
                                                            dfn$ms1 |
                                                            bin$mzmed %in% dfn$ms2,])
                        result <- result[rownames(unique(result[,
                                                                c("mzmed", "rtmed")])),]
                }
        }
        result$sm <- result$mzmed * massdiff
        result$smd <- ceiling(result$sm) - result$sm
        return(result)
}
#' Backup objects for furthor analysis
#' @param xset xcmsSet object
getbackup <- function(xset){

}
#' Function for Good/Bad peaks balance
#' @param xset xcmsSet object
#' @param cutoff minium point numbers for good peaks
#' @return Good peaks proportion and bad peaks proportion
#' @export
getprop <- function(xset, cutoff = 10){
        rt <- mean(diff(c(xset@rt$corrected[[1]])))
        peakwidthcutoff <- cutoff*rt
        t1 <- xcms::peaks(xset)
        index <- Reduce(c,xset@groupidx)
        t1s <- t1[index,]
        peakdiff <- t1s[,6]-t1s[,5]
        peakdiff0 <- t1[,6]-t1[,5]
        p <- sum(peakdiff>peakwidthcutoff)/nrow(t1)
        q <- sum(peakdiff<peakwidthcutoff)/sum(peakdiff0<peakwidthcutoff)
        return(c(p,q))
}
