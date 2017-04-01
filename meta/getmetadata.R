library(xcms)
library(BiocParallel)

getdata <-
        function(path,
                 index = F,
                 BPPARAM = SnowParam(workers = 12),
                 pmethod = 'hplcorbitrap',
                 ...) {
                cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
                if (index) {
                        cdffiles <- cdffiles[index]
                }
                if (pmethod == 'hplcorbitrap') {
                        xset <-
                                xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 2.5,
                                        peakwidth = c(10, 60),
                                        prefilter = c(3, 5000),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- group(xset,
                                              bw = 5,
                                              mzwid = 0.015)
                                xset2 <- retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        group(xset2,
                                              bw = 5,
                                              mzwid = 0.015)
                                xset3 <-
                                        fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'uplcorbitrap') {
                        xset <-
                                xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 2.5,
                                        peakwidth = c(5, 20),
                                        prefilter = c(3, 5000),
                                        ...
                                )
                        xset <- group(xset, bw = 2, mzwid = 0.015)
                        xset2 <- retcor(xset)
                        # you need group the peaks again for this corrected data
                        xset2 <- group(xset2, bw = 2, mzwid = 0.015)
                        xset3 <- fillPeaks(xset2, BPPARAM = BPPARAM)
                } else if (pmethod == 'hplcqtof') {
                        xset <-
                                xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 30,
                                        peakwidth = c(10, 60),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- group(xset,
                                              bw = 5,
                                              mzwid = 0.025)
                                xset2 <- retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        group(xset2,
                                              bw = 5,
                                              mzwid = 0.025)
                                xset3 <-
                                        fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'hplchqtof') {
                        xset <-
                                xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 15,
                                        peakwidth = c(10, 60),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- group(xset,
                                              bw = 5,
                                              mzwid = 0.015)
                                xset2 <- retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        group(xset2,
                                              bw = 5,
                                              mzwid = 0.015)
                                xset3 <-
                                        fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'uplcqtof') {
                        xset <-
                                xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 30,
                                        peakwidth = c(5, 20),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- group(xset,
                                              bw = 2,
                                              mzwid = 0.025)
                                xset2 <- retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        group(xset2,
                                              bw = 2,
                                              mzwid = 0.025)
                                xset3 <-
                                        fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'uplchqtof') {
                        xset <-
                                xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 15,
                                        peakwidth = c(5, 20),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- group(xset,
                                              bw = 2,
                                              mzwid = 0.015)
                                xset2 <- retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        group(xset2,
                                              bw = 2,
                                              mzwid = 0.015)
                                xset3 <-
                                        fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else{
                        xset <- xcmsSet(cdffiles, BPPARAM = BPPARAM, ...)
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- group(xset)
                                xset2 <- retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <- group(xset2)
                                xset3 <-
                                        fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                }
                return(xset3)
        }

getopmsdata <-
        function(path,
                 xsmethod = "matchedFilter",
                 fwhm = 35,
                 snthresh = 3,
                 step = 0.115,
                 steps = 3,
                 sigma = 14.8632580261593,
                 max = 5,
                 mzdiff = 0.455,
                 index = FALSE,
                 BPPARAM = SnowParam(workers = 12),
                 gmethod = "density",
                 bw = 12.4,
                 mzwid = 0.047,
                 minfrac = 0.892,
                 minsamp = 1,
                 gmax = 50,
                 rmethod = "obiwarp",
                 plottype = "none",
                 distFunc = "cor_opt",
                 profStep = 1,
                 center = 2,
                 response = 1,
                 gapInit = 0.4,
                 gapExtend = 2.064,
                 factorDiag = 2,
                 factorGap = 1,
                 localAlignment = 0,
                 ...) {
                cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
                xset <- xcmsSet(
                        cdffiles,
                        method = xsmethod,
                        fwhm = fwhm,
                        snthresh = snthresh,
                        step = step,
                        steps = steps,
                        sigma = sigma,
                        max = max,
                        mzdiff = mzdiff,
                        index = index,
                        BPPARAM = BPPARAM,
                        ...
                )
                xset <- group(xset)
                xset2 <- retcor(
                        xset,
                        method = rmethod,
                        plottype = plottype,
                        distFunc = distFunc,
                        profStep = profStep,
                        center = center,
                        response = response,
                        gapInit = gapInit,
                        gapExtend = gapExtend,
                        factorDiag = factorDiag,
                        factorGap = factorGap,
                        localAlignment = localAlignment
                )
                # you need group the peaks again for this corrected data
                xset2 <- group(
                        xset2,
                        method = gmethod,
                        bw = bw,
                        mzwid = mzwid,
                        minfrac = minfrac,
                        minsamp = minsamp,
                        max = gmax
                )
                xset3 <- fillPeaks(xset2, BPPARAM = BPPARAM)
                return(xset3)
        }

getopqedata <- function(path,
                        index = F,
                        xsmethod = "centWave",
                        peakwidth = c(14, 25),
                        ppm = 2.5,
                        noise = 0,
                        snthresh = 10,
                        mzdiff = -0.00395,
                        prefilter = c(3, 100),
                        mzCenterFun = "wMean",
                        integrate = 1,
                        fitgauss = FALSE,
                        verbose.columns = FALSE,
                        BPPARAM = SnowParam(workers = 12),
                        rmethod = "obiwarp",
                        plottype = "none",
                        distFunc = "cor_opt",
                        profStep = 1,
                        center = 2,
                        response = 1,
                        gapInit = 0.6176,
                        gapExtend = 2.4,
                        factorDiag = 2,
                        factorGap = 1,
                        localAlignment = 0,
                        gmethod = "density",
                        bw = 0.25,
                        mzwid = 0.0021748,
                        minfrac = 1,
                        minsamp = 1,
                        gmax = 50,
                        ...) {
        cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
        
        if (index) {
                cdffiles <- cdffiles[index]
        }
        xset <- xcmsSet(
                cdffiles,
                method = xsmethod,
                snthresh = snthresh,
                mzdiff = mzdiff,
                BPPARAM = BPPARAM,
                peakwidth = peakwidth,
                ppm = ppm,
                noise = noise,
                prefilter = prefilter,
                mzCenterFun = mzCenterFun,
                integrate = integrate,
                fitgauss = fitgauss,
                verbose.columns = verbose.columns,
                ...
        )
        if (index & length(index) == 1) {
                xset3 <- xset
        } else{
                xset <- group(
                        xset,
                        method = gmethod,
                        bw = bw,
                        mzwid = mzwid,
                        minfrac = minfrac,
                        minsamp = minsamp,
                        max = gmax
                )
                xset2 <- retcor(
                        xset,
                        method = rmethod,
                        plottype = plottype,
                        distFunc = distFunc,
                        profStep = profStep,
                        center = center,
                        response = response,
                        gapInit = gapInit,
                        gapExtend = gapExtend,
                        factorDiag = factorDiag,
                        factorGap = factorGap,
                        localAlignment = localAlignment
                )
                # you need group the peaks again for this corrected data
                xset2 <- group(
                        xset2,
                        method = gmethod,
                        bw = bw,
                        mzwid = mzwid,
                        minfrac = minfrac,
                        minsamp = minsamp,
                        max = gmax
                )
                xset3 <- fillPeaks(xset2, BPPARAM = BPPARAM)
        }
        return(xset3)
}

getupload <-
        function(xset,
                 method = "medret",
                 intensity = 'inio',
                 name = 'Peaklist') {
                peakIntensities <- groupval(xset, method, intensity)
                if (intensity == "intb") {
                        peakIntensities[is.na(peakIntensities)] = 0
                }
                data <-
                        rbind(group = as.character(phenoData(xset)$class), peakIntensities)
                filename <- paste0(name, '.csv')
                write.csv(data, file = filename)
                return(data)
        }

gettechrep <- function(xset,
                         method =  'medret',
                         intensity = 'into', ...){
        data <- t(groupval(xset, method, intensity, ...))
        lv <- xset@phenoData[, 1]
        mean <- aggregate(data, list(lv), mean)
        sd <- aggregate(data, list(lv), sd)
        suppressWarnings(rsd <- sd / mean * 100)
        result <-
                data.frame(cbind(t(mean[, -1]), t(sd[, -1]), t(rsd[, -1])))
        name <- unique(lv)
        colnames(result) <-
                c(paste0(name, 'mean'),
                  paste0(name, 'sd'),
                  paste0(name, 'rsd%'))
        datap <- groups(xset)
        report <- cbind.data.frame(datap, result)
        return(report)
}

gettechbiorep <-
        function(xset,
                 anno = F,
                 peaklist = F,
                 file = NULL,
                 method =  'medret',
                 intensity = 'into',
                 ...) {
                data <- t(groupval(xset, method, intensity, ...))
                lv <- xset@phenoData[, 1]
                lv2 <- xset@phenoData[, 2]
                mean <- aggregate(data, list(lv, lv2), mean)
                sd <- aggregate(data, list(lv, lv2), sd)
                suppressWarnings(rsd <- sd / mean * 100)
                result <-
                        data.frame(cbind(t(mean[, -c(1:2)]), t(sd[, -c(1:2)]), t(rsd[, -c(1:2)])))
                name <- unique(c(paste0(lv, lv2)))
                colnames(result) <-
                        c(paste0(name, 'mean'),
                          paste0(name, 'sd'),
                          paste0(name, 'rsd%'))
                datap <- groups(xset)
                report <- cbind.data.frame(datap, result)
                if (anno) {
                        anno <-
                                as.data.frame(cbind(xset@groups[, 1], xset@groups[, 4], cbind(t(
                                        mean[, -c(1:2)]
                                ))))
                        colnames(anno) <- c('mz', 'time', name)
                        return(anno)
                } else if (peaklist) {
                        result <- data.frame(t(mean[, -c(1:2)]))
                        data <- rbind(group = name, result)
                        write.csv(data, file = file)
                } else{
                        return(report)
                }
        }

getQCraw <- function(path, mzrange, rtrange, index = NULL) {
        cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
        if (index) {
                cdffiles <- cdffiles[index]
        }
        nsamples <- length(cdffiles)
        area <- numeric()
        for (i in 1:nsamples) {
                RAW <- xcmsRaw(cdffiles[i])
                peak <- rawEIC(RAW, mzrange, rtrange)
                area[i] <- sum(peak$intensity)
        }
        return(area)
}
