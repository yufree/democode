library(multtest)
library(xcms)
library(CAMERA)

getdata <- function(path,nSlaves=12,pmethod='hplcorbitrap',...){
  cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
  if(pmethod=='hplcorbitrap'){
    xset <- xcmsSet(cdffiles,nSlaves=nSlaves,method = "centWave",ppm=2.5,peakwidth=c(10,60),prefilter=c(3,5000),...)
    xset <- group(xset,bw=5,mzwid=0.015)
    xset2 <- retcor(xset)
    # you need group the peaks again for this corrected data
    xset2 <- group(xset2,bw=5,mzwid=0.015)
    xset3 <- fillPeaks(xset2,nSlaves=nSlaves)
  }else if(pmethod=='uplcorbitrap'){
    xset <- xcmsSet(cdffiles,nSlaves=nSlaves,method = "centWave",ppm=2.5,peakwidth=c(5,20),prefilter=c(3,5000),...)
    xset <- group(xset,bw=2,mzwid=0.015)
    xset2 <- retcor(xset)
    # you need group the peaks again for this corrected data
    xset2 <- group(xset2,bw=2,mzwid=0.015)
    xset3 <- fillPeaks(xset2,nSlaves=nSlaves)
  }else if(pmethod=='hplcqtof'){
    xset <- xcmsSet(cdffiles,nSlaves=nSlaves,method = "centWave",ppm=30,peakwidth=c(10,60),prefilter=c(0,0),...)
    xset <- group(xset,bw=5,mzwid=0.025)
    xset2 <- retcor(xset)
    # you need group the peaks again for this corrected data
    xset2 <- group(xset2,bw=5,mzwid=0.025)
    xset3 <- fillPeaks(xset2,nSlaves=nSlaves)
  }else if(pmethod=='hplchqtof'){
    xset <- xcmsSet(cdffiles,nSlaves=nSlaves,method = "centWave",ppm=15,peakwidth=c(10,60),prefilter=c(0,0),...)
    xset <- group(xset,bw=5,mzwid=0.015)
    xset2 <- retcor(xset)
    # you need group the peaks again for this corrected data
    xset2 <- group(xset2,bw=5,mzwid=0.015)
    xset3 <- fillPeaks(xset2,nSlaves=nSlaves)
  }else if(pmethod=='uplcqtof'){
    xset <- xcmsSet(cdffiles,nSlaves=nSlaves,method = "centWave",ppm=30,peakwidth=c(5,20),prefilter=c(0,0),...)
    xset <- group(xset,bw=2,mzwid=0.025)
    xset2 <- retcor(xset)
    # you need group the peaks again for this corrected data
    xset2 <- group(xset2,bw=2,mzwid=0.025)
    xset3 <- fillPeaks(xset2,nSlaves=nSlaves)
  }else if(pmethod=='uplchqtof'){
    xset <- xcmsSet(cdffiles,nSlaves=nSlaves,method = "centWave",ppm=15,peakwidth=c(5,20),prefilter=c(0,0),...)
    xset <- group(xset,bw=2,mzwid=0.015)
    xset2 <- retcor(xset)
    # you need group the peaks again for this corrected data
    xset2 <- group(xset2,bw=2,mzwid=0.015)
    xset3 <- fillPeaks(xset2,nSlaves=nSlaves)
  }else{
    xset <- xcmsSet(cdffiles,nSlaves=nSlaves,...)
    xset <- group(xset)
    xset2 <- retcor(xset)
    # you need group the peaks again for this corrected data
    xset2 <- group(xset2)
    xset3 <- fillPeaks(xset2,nSlaves=nSlaves)
  }
  return(xset3)
}

getopmsdata <- function(path,xsmethod = "matchedFilter",fwhm=35,snthresh=3, step=0.115, steps=3, sigma=14.8632580261593, max=5, mzdiff=0.455, index=FALSE, nSlaves=12,gmethod="density", 
                    bw=12.4, mzwid=0.047, minfrac=0.892, minsamp=1, gmax=50,rmethod="obiwarp",
                    plottype="none", distFunc="cor_opt", profStep=1, center=2, response=1, gapInit=0.4, gapExtend=2.064, factorDiag=2, factorGap=1, localAlignment=0,...){
  cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
  xset <- xcmsSet(cdffiles,
                  method=xsmethod, 
                  fwhm=fwhm, 
                  snthresh=snthresh, 
                  step=step, 
                  steps=steps, 
                  sigma=sigma, 
                  max=max, 
                  mzdiff=mzdiff, 
                  index=index, 
                  nSlaves=nSlaves,...)
  xset <- group(xset)
  xset2 <- retcor(xset, method=rmethod,
                  plottype=plottype,
                  distFunc=distFunc, 
                  profStep=profStep, 
                  center=center, 
                  response=response,
                  gapInit=gapInit, 
                  gapExtend=gapExtend, 
                  factorDiag=factorDiag, 
                  factorGap=factorGap, 
                  localAlignment=localAlignment)
  # you need group the peaks again for this corrected data
  xset2 <- group(xset2,
                 method=gmethod, 
                 bw=bw,
                 mzwid=mzwid, 
                 minfrac=minfrac, 
                 minsamp=minsamp, 
                 max=gmax)
  xset3 <- fillPeaks(xset2,nSlaves=nSlaves)
  return(xset3)
}

getopqedata <- function(path,xsmethod = "centWave", 
                        peakwidth=c(14, 25), ppm=2.5, 
                        noise=0, snthresh=10, 
                        mzdiff=-0.00395, prefilter=c(3, 100),
                        mzCenterFun="wMean", integrate=1,
                        fitgauss=FALSE, verbose.columns=FALSE,
                        nSlaves=12,rmethod="obiwarp",
                        plottype="none", distFunc="cor_opt", 
                        profStep=1, center=2, 
                        response=1, gapInit=0.6176, 
                        gapExtend=2.4, factorDiag=2, 
                        factorGap=1, localAlignment=0,
                        gmethod="density", bw=0.25,
                        mzwid=0.0021748, minfrac=1,
                        minsamp=1, gmax=50,...){
  cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
  xset <- xcmsSet(cdffiles,
                  method=xsmethod, 
                  snthresh=snthresh, 
                  mzdiff=mzdiff,
                  nSlaves=nSlaves,
                  peakwidth = peakwidth,
                  ppm=ppm, 
                  noise=noise, 
                  prefilter=prefilter, 
                  mzCenterFun=mzCenterFun, 
                  integrate=integrate, 
                  fitgauss=fitgauss, 
                  verbose.columns=verbose.columns,
                  ...)
  xset <- group(xset)
  xset2 <- retcor(xset, method=rmethod,
                  plottype=plottype,
                  distFunc=distFunc, 
                  profStep=profStep, 
                  center=center, 
                  response=response,
                  gapInit=gapInit, 
                  gapExtend=gapExtend, 
                  factorDiag=factorDiag, 
                  factorGap=factorGap, 
                  localAlignment=localAlignment)
  # you need group the peaks again for this corrected data
  xset2 <- group(xset2,
                 method=gmethod, 
                 bw=bw,
                 mzwid=mzwid, 
                 minfrac=minfrac, 
                 minsamp=minsamp, 
                 max=gmax)
  xset3 <- fillPeaks(xset2,nSlaves=nSlaves)
  return(xset3)
}
