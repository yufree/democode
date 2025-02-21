## ----setup, include=FALSE----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----------------------------------------------------------------------------------------------------
library(mzrtsim)
# load the high resolution MS1 database from MoNA
data(monahrms1)
# Set peak height
ph1 <- c(rep(5,30),rep(10,40),rep(15,30))
ph2 <- c(rep(5,20),rep(10,30),rep(15,50))
# Set retention time
rt <- seq(10,590,length.out=100)
# rt[c(21:30,51:70)]
# display peaks profile
plot(ph1~rt,cex=0.5,col='blue')
points(ph2~rt,col='red')
# for reproducible purpose
set.seed(1)
# select compounds
compound <- sample(c(1:1114),100)
set.seed(2)
# sample response factors for each compounds
rf <- sample(c(100:10000),100)

# 30% changed
for(i in c(1:10)){
  # with unique, only one spectra will be used for simulated compound and no duplicated spectra for the same compound as database will contain multiple spectra for the same compound
  simmzml(name=paste0('sim/case/case',i),db=monahrms1,pheight = ph1,compound=compound,rtime = rt, rf=rf,unique = T,matrix = T)
}

for(i in c(1:10)){
  # set different peak width
  simmzml(name=paste0('sim/control/control',i),db=monahrms1,pheight = ph2,compound=compound,rtime = rt, rf=rf,unique = T,matrix = T)
}
files <- list.files('sim',full.names = T,recursive = T,pattern = '*.csv')
name <- basename(files)
file.rename(files,paste0('simcsv/',name))


## ----------------------------------------------------------------------------------------------------
library(mzrtsim)
# the following part is the same with previous simulation
data(monahrms1)
ph1 <- c(rep(5,30),rep(10,40),rep(15,30))
ph2 <- c(rep(5,20),rep(10,30),rep(15,50))
rt <- seq(10,590,length.out=100)
rt[c(21:30,51:70)]
plot(ph1~rt,cex=0.5,col='blue')
points(ph2~rt,col='red')
set.seed(1)
compound <- sample(c(1:1114),100)
set.seed(2)
SNR <- sample(c(100:10000),100)

# 30% changed
for(i in c(1:10)){
  # set intenisty cutoff as 0
  simmzml(name=paste0('sim2/case/case',i),db=monahrms1,pheight = ph1,compound=compound,rtime = rt, rf=rf,unique = T,inscutoff = 0,matrix = T)
}

for(i in c(1:10)){
  # set intenisty cutoff as 0
  simmzml(name=paste0('sim2/control/control',i),db=monahrms1,pheight = ph2,compound=compound,rtime = rt, rf=rf,unique = T,inscutoff = 0,matrix = T)
}

files <- list.files('sim2',full.names = T,recursive = T,pattern = '*.csv')
name <- basename(files)
file.rename(files,paste0('simcsv2/',name))


## ----------------------------------------------------------------------------------------------------
library(mzrtsim)
# the following part is the same with previous simulation
data(monahrms1)
ph <- c(rep(5,30),rep(10,40),rep(15,30))
rt <- seq(10,590,length.out=100)
plot(ph1~rt,cex=0.5,col='blue')
points(ph2~rt,col='red')
set.seed(1)
compound <- sample(c(1:1114),100)
set.seed(2)
SNR <- sample(c(100:10000),100)

for(i in c(1:10)){
  # change tailing factor to 1.5 to simulate tailing peaks, add matrix effect
  simmzml(name=paste0('sim3/tailing/tailing',i),db=monahrms1,pheight = ph,compound=compound,rtime = rt, SNR=SNR,unique = T,tailingfactor = 1.5,matrix = T)
}

for(i in c(1:10)){
  # change tailing factor to 1 to simulate normal peaks, add matrix effect
  simmzml(name=paste0('sim3/normal/normal',i),db=monahrms1,pheight = ph,compound=compound,rtime = rt, SNR=SNR,unique = T,tailingfactor = 1,matrix = T)
}

for(i in c(1:10)){
  # change tailing factor to 0.8 to simulate leading peaks, add matrix effect
  simmzml(name=paste0('sim3/leading/leading',i),db=monahrms1,pwidth = ph,compound=compound,rtime = rt, SNR=SNR,unique = T,tailingfactor = 0.8,matrix = T)
}

files <- list.files('sim3',full.names = T,recursive = T,pattern = '*.csv')
name <- basename(files)
file.rename(files,paste0('simcsv3/',name))

## ----mzmine------------------------------------------------------------------------------------------
# GUI mzMine 4.5 is used to extract peaks. Raw data is imported as mzML files and then mass detection with default setting is performed. Feature detection used chromatogram builder in LC-MS with default setting. When feature detection is done for each file, alignment was perfomred with join aligner with weights of m/z 0.5 and retention time 0.5. The aligned feature list is exported for further comparison.

## ----------------------------------------------------------------------------------------------------
library(xcms)
library(MsExperiment)
# peak picking for simulated peaks
path <- 'sim'
files <- list.files(path,pattern = ".CDF|.mzXML|.mzML",full.names = T,recursive = T)
group <- xcms:::phenoDataFromPaths(files)
  if(NCOL(group)==1){
      sample_group <- group$class
  }else{
      cols <- colnames(group)
      sample_group <-  do.call(paste,c(group[cols],sep='_'))
  }
  sample_name=sub(basename(files),pattern = ".CDF|.mzXML|.mzML",replacement = '')
  pd <- cbind.data.frame(sample_name, sample_group)
  
raw_data <- readMsExperiment(spectraFiles = files, sampleData = pd)
# change peak width to [5,15] to cover simulated compound
cwp <- CentWaveParam(ppm=5,peakwidth = c(5,15))
xdata <- findChromPeaks(raw_data, param = cwp)
xdata <- groupChromPeaks(xdata, param = PeakDensityParam(pd$sample_group))
xdata<- fillChromPeaks(xdata, param = ChromPeakAreaParam())
# extract peaks profile as csv file
re <- featureDefinitions(xdata)[, c("mzmed","rtmed")]
# 340 peaks
# save peaks profile with name
compound <- read.csv('simcsv/case1.csv')
# 593 peaks
# simulate mass range [100,1000] for MS1 scan
compoundsub <- compound[compound$mz>100&compound$mz<1000,]
# 533 peaks
# align simulated peaks with detected peaks
align2 <- enviGCMS::getalign(compoundsub$mz,re$mzmed,compoundsub$rt,re$rtmed)
length(unique(compound$name[align2$xid]))
# 74
length(unique(paste(align2$mz2,align2$rt2)))
# 330 found peaks match to 335 peaks
library(SummarizedExperiment)
res <- quantify(xdata, value = "maxo", method = "max")
data <- assay(res)
library(genefilter)
# check changed peaks
rda <- rowttests(as.matrix(data),fac=as.factor(pd$sample_group))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05,na.rm = T)
# 134
df <- cbind.data.frame(re,data)
write.csv(df,'simxcms.csv')


## ----------------------------------------------------------------------------------------------------
library(xcms)
library(MsExperiment)
# add this for mac os
# setMSnbaseFastLoad(opt = F)
# peak picking for simulated peaks
path <- 'sim2'
files <- list.files(path,pattern = ".CDF|.mzXML|.mzML",full.names = T,recursive = T)
group <- xcms:::phenoDataFromPaths(files)
  if(NCOL(group)==1){
      sample_group <- group$class
  }else{
      cols <- colnames(group)
      sample_group <-  do.call(paste,c(group[cols],sep='_'))
  }
  sample_name=sub(basename(files),pattern = ".CDF|.mzXML|.mzML",replacement = '')
  pd <- cbind.data.frame(sample_name, sample_group)
  
raw_data <- readMsExperiment(spectraFiles = files, sampleData = pd)
# change peak width to [5,15] to cover simulated compound
cwp <- CentWaveParam(ppm=5,peakwidth = c(5,15))
xdata <- findChromPeaks(raw_data, param = cwp)
xdata <- groupChromPeaks(xdata, param = PeakDensityParam(pd$sample_group))
xdata<- fillChromPeaks(xdata, param = ChromPeakAreaParam())
# extract peaks profile as csv file
re <- featureDefinitions(xdata)[, c("mzmed","rtmed")]
# 2099 peaks
# save peaks profile with name
compound <- read.csv('simcsv2/case1.csv')
# 593 peaks
# simulate mass range [100,1000] for MS1 scan
compoundsub <- compound[compound$mz>100&compound$mz<1000,]
# 533 peaks
# align simulated peaks with detected peaks
align2 <- enviGCMS::getalign(compoundsub$mz,re$mzmed,compoundsub$rt,re$rtmed)
length(unique(compound$name[align2$xid]))
# 65
length(unique(paste(align2$mz2,align2$rt2)))
# 1972 found peaks match to 2051 peaks
library(SummarizedExperiment)
res <- quantify(xdata, value = "maxo", method = "max")
data <- assay(res)
library(genefilter)
# check changed peaks
rda <- rowttests(as.matrix(data),fac=as.factor(pd$sample_group))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05,na.rm = T)
# 916
df <- cbind.data.frame(re,data)
write.csv(df,'sim2xcms.csv')


## ----leading-----------------------------------------------------------------------------------------
library(xcms)
library(MsExperiment)
# add this for mac os
# setMSnbaseFastLoad(opt = F)
# peak picking for simulated peaks
path <- 'sim3/leading/'
files <- list.files(path,pattern = ".CDF|.mzXML|.mzML",full.names = T,recursive = T)
group <- xcms:::phenoDataFromPaths(files)
  if(NCOL(group)==1){
      sample_group <- group$class
  }else{
      cols <- colnames(group)
      sample_group <-  do.call(paste,c(group[cols],sep='_'))
  }
  sample_name=sub(basename(files),pattern = ".CDF|.mzXML|.mzML",replacement = '')
  pd <- cbind.data.frame(sample_name, sample_group)
  
raw_data <- readMsExperiment(spectraFiles = files, sampleData = pd)
# change peak width to [5,15] to cover simulated compound
cwp <- CentWaveParam(ppm=5,peakwidth = c(5,15))
xdata <- findChromPeaks(raw_data, param = cwp)
xdata <- groupChromPeaks(xdata, param = PeakDensityParam(pd$sample_group))
xdata<- fillChromPeaks(xdata, param = ChromPeakAreaParam())
# extract peaks profile as csv file
re <- featureDefinitions(xdata)[, c("mzmed","rtmed")]
# 412 peaks
# save peaks profile with name
compound <- read.csv('simcsv3/leading1.csv')
# 593 peaks
# simulate mass range [100,1000] for MS1 scan
compoundsub <- compound[compound$mz>100&compound$mz<1000,]
# 533 peaks
# align simulated peaks with detected peaks
align2 <- enviGCMS::getalign(compoundsub$mz,re$mzmed,compoundsub$rt,re$rtmed,ppm = 5,deltart = 5)
length(unique(compound$name[align2$xid]))
# 83
length(unique(paste(align2$mz2,align2$rt2)))
# 396 found peaks match to 397 peaks
library(SummarizedExperiment)
res <- quantify(xdata, value = "maxo", method = "max")
data <- assay(res)
leading <- cbind.data.frame(re,data)
write.csv(leading,'sim3xcmsleading.csv')


## ----normal------------------------------------------------------------------------------------------
library(xcms)
library(MsExperiment)
# add this for mac os
# setMSnbaseFastLoad(opt = F)
# peak picking for simulated peaks
path <- 'sim3/normal'
files <- list.files(path,pattern = ".CDF|.mzXML|.mzML",full.names = T,recursive = T)
group <- xcms:::phenoDataFromPaths(files)
  if(NCOL(group)==1){
      sample_group <- group$class
  }else{
      cols <- colnames(group)
      sample_group <-  do.call(paste,c(group[cols],sep='_'))
  }
  sample_name=sub(basename(files),pattern = ".CDF|.mzXML|.mzML",replacement = '')
  pd <- cbind.data.frame(sample_name, sample_group)
  
raw_data <- readMsExperiment(spectraFiles = files, sampleData = pd)
# change peak width to [5,15] to cover simulated compound
cwp <- CentWaveParam(ppm=5,peakwidth = c(5,15))
xdata <- findChromPeaks(raw_data, param = cwp)
xdata <- groupChromPeaks(xdata, param = PeakDensityParam(pd$sample_group))
xdata<- fillChromPeaks(xdata, param = ChromPeakAreaParam())
# extract peaks profile as csv file
re <- featureDefinitions(xdata)[, c("mzmed","rtmed")]
# 377 peaks
# save peaks profile with name
compound <- read.csv('simcsv3/normal1.csv')
# 593 peaks
# simulate mass range [100,1000] for MS1 scan
compoundsub <- compound[compound$mz>100&compound$mz<1000,]
# 533 peaks
# align simulated peaks with detected peaks
align2 <- enviGCMS::getalign(compoundsub$mz,re$mzmed,compoundsub$rt,re$rtmed,ppm = 5,deltart = 5)
length(unique(compoundsub$name[align2$xid]))
# 82
length(unique(paste(align2$mz2,align2$rt2)))
# 367 found peaks match to 367 peaks
library(SummarizedExperiment)
res <- quantify(xdata, value = "maxo", method = "max")
data <- assay(res)
normal <- cbind.data.frame(re,data)
write.csv(normal,'sim3xcmsnormal.csv')


## ----tailing-----------------------------------------------------------------------------------------
library(xcms)
library(MsExperiment)
# add this for mac os
# setMSnbaseFastLoad(opt = F)
# peak picking for simulated peaks
path <- 'sim3/tailing'
files <- list.files(path,pattern = ".CDF|.mzXML|.mzML",full.names = T,recursive = T)
group <- xcms:::phenoDataFromPaths(files)
  if(NCOL(group)==1){
      sample_group <- group$class
  }else{
      cols <- colnames(group)
      sample_group <-  do.call(paste,c(group[cols],sep='_'))
  }
  sample_name=sub(basename(files),pattern = ".CDF|.mzXML|.mzML",replacement = '')
  pd <- cbind.data.frame(sample_name, sample_group)
  
raw_data <- readMsExperiment(spectraFiles = files, sampleData = pd)
# change peak width to [5,15] to cover simulated compound
cwp <- CentWaveParam(ppm=5,peakwidth = c(5,15))
xdata <- findChromPeaks(raw_data, param = cwp)
xdata <- groupChromPeaks(xdata, param = PeakDensityParam(pd$sample_group))
xdata<- fillChromPeaks(xdata, param = ChromPeakAreaParam())
# extract peaks profile as csv file
re <- featureDefinitions(xdata)[, c("mzmed","rtmed")]
# 294 peaks
# save peaks profile with name
compound <- read.csv('simcsv3/tailing1.csv')
# 593 peaks
# simulate mass range [100,1000] for MS1 scan
compoundsub <- compound[compound$mz>100&compound$mz<1000,]
# 533 peaks
# align simulated peaks with detected peaks
align2 <- enviGCMS::getalign(compoundsub$mz,re$mzmed,compoundsub$rt,re$rtmed,ppm = 5,deltart = 5)
length(unique(compound$name[align2$xid]))
# 77
length(unique(paste(align2$mz2,align2$rt2)))
# 289 found peaks match to 289 peaks
library(SummarizedExperiment)
res <- quantify(xdata, value = "maxo", method = "max")
data <- assay(res)
tailing <- cbind.data.frame(re,data)
write.csv(tailing,'sim3xcmstailing.csv')


## ----------------------------------------------------------------------------------------------------
# you might need to install python and pyopenms to run the following code
reticulate::use_python('/opt/homebrew/Caskroom/miniconda/base/bin/python')


## import os
## import shutil
## import requests
## import pandas as pd
## from pyopenms import *
## import os
## import numpy as np
## file1 = "sim/case/"
## file2 = "sim/control/"
## 
## mzML_files = [file1+x for x in os.listdir(file1)]+[file2+x for x in os.listdir(file2)]
## 
## feature_maps = []
## for file in mzML_files:
##     # load mzML file into MSExperiment
##     exp = MSExperiment()
##     MzMLFile().load(file, exp) # load each mzML file to an OpenMS file format (MSExperiment)
## 
##     # mass trace detection
##     mass_traces = [] # introduce an empty list where the mass traces will be loaded
##     mtd = MassTraceDetection()
##     mtd_par = mtd.getDefaults() # get the default parameters in order to edit them
##     mtd_par.setValue("mass_error_ppm", 5.0) # high-res instrument, orbitraps
##     mtd_par.setValue("noise_threshold_int", 1.0e02) # data-dependent (usually works for orbitraps)
##     mtd.setParameters(mtd_par) # set the new parameters
##     mtd.run(exp, mass_traces, 0) # run mass trace detection
## 
##     # elution peak detection
##     mass_traces_deconvol = []
##     epd = ElutionPeakDetection()
##     epd_par = epd.getDefaults()
##     epd_par.setValue("width_filtering", "fixed") # The fixed setting filters out mass traces outside the [min_fwhm: 1.0, max_fwhm: 60.0] interval
##     epd.setParameters(epd_par)
##     epd.detectPeaks(mass_traces, mass_traces_deconvol)
## 
##     # feature detection
##     feature_map = FeatureMap() # output features
##     chrom_out = [] # output chromatograms
##     ffm = FeatureFindingMetabo()
##     ffm_par = ffm.getDefaults()
##     #ffm_par.setValue("remove_single_traces", "true") # remove mass traces without satellite isotopic traces
##     #ffm.setParameters(ffm_par)
##     ffm.run(mass_traces_deconvol, feature_map, chrom_out)
##     feature_map.setUniqueIds() # Assigns a new, valid unique id per feature
##     feature_map.setPrimaryMSRunPath([file.encode()]) # Sets the file path to the primary MS run (usually the mzML file)
##     feature_maps.append(feature_map)
## 
## 
## # use as reference for alignment, the file with the largest number of features (works well if you have a pooled QC for example)
## 
## ref_index = feature_maps.index(sorted(feature_maps, key=lambda x: x.size())[-1])
## 
## aligner = MapAlignmentAlgorithmPoseClustering()
## 
## trafos = {}
## 
## # parameter optimization
## aligner_par= aligner.getDefaults()
## aligner_par.setValue("max_num_peaks_considered", -1) # infinite
## aligner_par.setValue("pairfinder:distance_MZ:max_difference", 10.0) # Never pair features with larger m/z distance
## aligner_par.setValue("pairfinder:distance_MZ:unit", "ppm")
## aligner.setParameters(aligner_par)
## aligner.setReference(feature_maps[ref_index])
## 
## for feature_map in feature_maps[:ref_index] + feature_maps[ref_index+1:]:
##     trafo = TransformationDescription() # save the transformed data points
##     aligner.align(feature_map, trafo)
##     trafos[feature_map.getMetaValue("spectra_data")[0].decode()] = trafo
##     transformer = MapAlignmentTransformer()
##     transformer.transformRetentionTimes(feature_map, trafo, True)
## 
## 
## feature_grouper = FeatureGroupingAlgorithmKD()
## 
## consensus_map = ConsensusMap()
## file_descriptions = consensus_map.getColumnHeaders()
## 
## for i, feature_map in enumerate(feature_maps):
##     file_description = file_descriptions.get(i, ColumnHeader())
##     file_description.filename = os.path.basename(
##         feature_map.getMetaValue("spectra_data")[0].decode())
##     file_description.size = feature_map.size()
##     file_descriptions[i] = file_description
## 
## 
## feature_grouper.group(feature_maps, consensus_map)
## consensus_map.setColumnHeaders(file_descriptions)
## consensus_map.setUniqueIds()
## ConsensusXMLFile().store("FeatureMatrix.consensusXML", consensus_map)
## df = consensus_map.get_df()
## df.to_csv('simopenms.csv')

## import os
## import shutil
## import requests
## import pandas as pd
## from pyopenms import *
## import os
## import numpy as np
## file1 = "sim2/case/"
## file2 = "sim2/control/"
## 
## mzML_files = [file1+x for x in os.listdir(file1)]+[file2+x for x in os.listdir(file2)]
## 
## feature_maps = []
## for file in mzML_files:
##     # load mzML file into MSExperiment
##     exp = MSExperiment()
##     MzMLFile().load(file, exp) # load each mzML file to an OpenMS file format (MSExperiment)
## 
##     # mass trace detection
##     mass_traces = [] # introduce an empty list where the mass traces will be loaded
##     mtd = MassTraceDetection()
##     mtd_par = mtd.getDefaults() # get the default parameters in order to edit them
##     mtd_par.setValue("mass_error_ppm", 5.0) # high-res instrument, orbitraps
##     mtd_par.setValue("noise_threshold_int", 1.0e02) # data-dependent (usually works for orbitraps)
##     mtd.setParameters(mtd_par) # set the new parameters
##     mtd.run(exp, mass_traces, 0) # run mass trace detection
## 
##     # elution peak detection
##     mass_traces_deconvol = []
##     epd = ElutionPeakDetection()
##     epd_par = epd.getDefaults()
##     epd_par.setValue("width_filtering", "fixed") # The fixed setting filters out mass traces outside the [min_fwhm: 1.0, max_fwhm: 60.0] interval
##     epd.setParameters(epd_par)
##     epd.detectPeaks(mass_traces, mass_traces_deconvol)
## 
##     # feature detection
##     feature_map = FeatureMap() # output features
##     chrom_out = [] # output chromatograms
##     ffm = FeatureFindingMetabo()
##     ffm_par = ffm.getDefaults()
##     #ffm_par.setValue("remove_single_traces", "true") # remove mass traces without satellite isotopic traces
##     #ffm.setParameters(ffm_par)
##     ffm.run(mass_traces_deconvol, feature_map, chrom_out)
##     feature_map.setUniqueIds() # Assigns a new, valid unique id per feature
##     feature_map.setPrimaryMSRunPath([file.encode()]) # Sets the file path to the primary MS run (usually the mzML file)
##     feature_maps.append(feature_map)
## 
## 
## # use as reference for alignment, the file with the largest number of features (works well if you have a pooled QC for example)
## 
## ref_index = feature_maps.index(sorted(feature_maps, key=lambda x: x.size())[-1])
## 
## aligner = MapAlignmentAlgorithmPoseClustering()
## 
## trafos = {}
## 
## # parameter optimization
## aligner_par= aligner.getDefaults()
## aligner_par.setValue("max_num_peaks_considered", -1) # infinite
## aligner_par.setValue("pairfinder:distance_MZ:max_difference", 10.0) # Never pair features with larger m/z distance
## aligner_par.setValue("pairfinder:distance_MZ:unit", "ppm")
## aligner.setParameters(aligner_par)
## aligner.setReference(feature_maps[ref_index])
## 
## for feature_map in feature_maps[:ref_index] + feature_maps[ref_index+1:]:
##     trafo = TransformationDescription() # save the transformed data points
##     aligner.align(feature_map, trafo)
##     trafos[feature_map.getMetaValue("spectra_data")[0].decode()] = trafo
##     transformer = MapAlignmentTransformer()
##     transformer.transformRetentionTimes(feature_map, trafo, True)
## 
## 
## feature_grouper = FeatureGroupingAlgorithmKD()
## 
## consensus_map = ConsensusMap()
## file_descriptions = consensus_map.getColumnHeaders()
## 
## for i, feature_map in enumerate(feature_maps):
##     file_description = file_descriptions.get(i, ColumnHeader())
##     file_description.filename = os.path.basename(
##         feature_map.getMetaValue("spectra_data")[0].decode())
##     file_description.size = feature_map.size()
##     file_descriptions[i] = file_description
## 
## 
## feature_grouper.group(feature_maps, consensus_map)
## consensus_map.setColumnHeaders(file_descriptions)
## consensus_map.setUniqueIds()
## ConsensusXMLFile().store("FeatureMatrix.consensusXML", consensus_map)
## df = consensus_map.get_df()
## df.to_csv('sim2openms.csv')

## import os
## import shutil
## import requests
## import pandas as pd
## from pyopenms import *
## import os
## import numpy as np
## file1 = "sim3/leading/"
## file2 = "sim3/normal/"
## file3 = "sim3/tailing/"
## 
## # mzML_files = [file1+x for x in os.listdir(file1)]+[file2+x for x in os.listdir(file2)]+[file3+x for x in os.listdir(file3)]
## 
## # mzML_files = [file2+x for x in os.listdir(file2)]+[file3+x for x in os.listdir(file3)]
## 
## # mzML_files = [file3+x for x in os.listdir(file3)]
## 
## mzML_files = [file1+x for x in os.listdir(file1)]
## 
## feature_maps = []
## for file in mzML_files:
##     # load mzML file into MSExperiment
##     exp = MSExperiment()
##     MzMLFile().load(file, exp) # load each mzML file to an OpenMS file format (MSExperiment)
## 
##     # mass trace detection
##     mass_traces = [] # introduce an empty list where the mass traces will be loaded
##     mtd = MassTraceDetection()
##     mtd_par = mtd.getDefaults() # get the default parameters in order to edit them
##     mtd_par.setValue("mass_error_ppm", 5.0) # high-res instrument, orbitraps
##     mtd_par.setValue("noise_threshold_int", 1.0e02) # data-dependent (usually works for orbitraps)
##     mtd.setParameters(mtd_par) # set the new parameters
##     mtd.run(exp, mass_traces, 0) # run mass trace detection
## 
##     # elution peak detection
##     mass_traces_deconvol = []
##     epd = ElutionPeakDetection()
##     epd_par = epd.getDefaults()
##     epd_par.setValue("width_filtering", "fixed") # The fixed setting filters out mass traces outside the [min_fwhm: 1.0, max_fwhm: 60.0] interval
##     epd.setParameters(epd_par)
##     epd.detectPeaks(mass_traces, mass_traces_deconvol)
## 
##     # feature detection
##     feature_map = FeatureMap() # output features
##     chrom_out = [] # output chromatograms
##     ffm = FeatureFindingMetabo()
##     ffm_par = ffm.getDefaults()
##     #ffm_par.setValue("remove_single_traces", "true") # remove mass traces without satellite isotopic traces
##     #ffm.setParameters(ffm_par)
##     ffm.run(mass_traces_deconvol, feature_map, chrom_out)
##     feature_map.setUniqueIds() # Assigns a new, valid unique id per feature
##     feature_map.setPrimaryMSRunPath([file.encode()]) # Sets the file path to the primary MS run (usually the mzML file)
##     feature_maps.append(feature_map)
## 
## 
## # use as reference for alignment, the file with the largest number of features (works well if you have a pooled QC for example)
## 
## ref_index = feature_maps.index(sorted(feature_maps, key=lambda x: x.size())[-1])
## 
## aligner = MapAlignmentAlgorithmPoseClustering()
## 
## trafos = {}
## 
## # parameter optimization
## aligner_par= aligner.getDefaults()
## aligner_par.setValue("max_num_peaks_considered", -1) # infinite
## aligner_par.setValue("pairfinder:distance_MZ:max_difference", 10.0) # Never pair features with larger m/z distance
## aligner_par.setValue("pairfinder:distance_MZ:unit", "ppm")
## aligner.setParameters(aligner_par)
## aligner.setReference(feature_maps[ref_index])
## 
## for feature_map in feature_maps[:ref_index] + feature_maps[ref_index+1:]:
##     trafo = TransformationDescription() # save the transformed data points
##     aligner.align(feature_map, trafo)
##     trafos[feature_map.getMetaValue("spectra_data")[0].decode()] = trafo
##     transformer = MapAlignmentTransformer()
##     transformer.transformRetentionTimes(feature_map, trafo, True)
## 
## 
## feature_grouper = FeatureGroupingAlgorithmKD()
## 
## consensus_map = ConsensusMap()
## file_descriptions = consensus_map.getColumnHeaders()
## 
## for i, feature_map in enumerate(feature_maps):
##     file_description = file_descriptions.get(i, ColumnHeader())
##     file_description.filename = os.path.basename(
##         feature_map.getMetaValue("spectra_data")[0].decode())
##     file_description.size = feature_map.size()
##     file_descriptions[i] = file_description
## 
## 
## feature_grouper.group(feature_maps, consensus_map)
## consensus_map.setColumnHeaders(file_descriptions)
## consensus_map.setUniqueIds()
## ConsensusXMLFile().store("FeatureMatrix.consensusXML", consensus_map)
## df = consensus_map.get_df()
## df.to_csv('sim3openmsleading.csv')
## 
## mzML_files = [file2+x for x in os.listdir(file2)]
## 
## feature_maps = []
## for file in mzML_files:
##     # load mzML file into MSExperiment
##     exp = MSExperiment()
##     MzMLFile().load(file, exp) # load each mzML file to an OpenMS file format (MSExperiment)
## 
##     # mass trace detection
##     mass_traces = [] # introduce an empty list where the mass traces will be loaded
##     mtd = MassTraceDetection()
##     mtd_par = mtd.getDefaults() # get the default parameters in order to edit them
##     mtd_par.setValue("mass_error_ppm", 5.0) # high-res instrument, orbitraps
##     mtd_par.setValue("noise_threshold_int", 1.0e02) # data-dependent (usually works for orbitraps)
##     mtd.setParameters(mtd_par) # set the new parameters
##     mtd.run(exp, mass_traces, 0) # run mass trace detection
## 
##     # elution peak detection
##     mass_traces_deconvol = []
##     epd = ElutionPeakDetection()
##     epd_par = epd.getDefaults()
##     epd_par.setValue("width_filtering", "fixed") # The fixed setting filters out mass traces outside the [min_fwhm: 1.0, max_fwhm: 60.0] interval
##     epd.setParameters(epd_par)
##     epd.detectPeaks(mass_traces, mass_traces_deconvol)
## 
##     # feature detection
##     feature_map = FeatureMap() # output features
##     chrom_out = [] # output chromatograms
##     ffm = FeatureFindingMetabo()
##     ffm_par = ffm.getDefaults()
##     #ffm_par.setValue("remove_single_traces", "true") # remove mass traces without satellite isotopic traces
##     #ffm.setParameters(ffm_par)
##     ffm.run(mass_traces_deconvol, feature_map, chrom_out)
##     feature_map.setUniqueIds() # Assigns a new, valid unique id per feature
##     feature_map.setPrimaryMSRunPath([file.encode()]) # Sets the file path to the primary MS run (usually the mzML file)
##     feature_maps.append(feature_map)
## 
## 
## # use as reference for alignment, the file with the largest number of features (works well if you have a pooled QC for example)
## 
## ref_index = feature_maps.index(sorted(feature_maps, key=lambda x: x.size())[-1])
## 
## aligner = MapAlignmentAlgorithmPoseClustering()
## 
## trafos = {}
## 
## # parameter optimization
## aligner_par= aligner.getDefaults()
## aligner_par.setValue("max_num_peaks_considered", -1) # infinite
## aligner_par.setValue("pairfinder:distance_MZ:max_difference", 10.0) # Never pair features with larger m/z distance
## aligner_par.setValue("pairfinder:distance_MZ:unit", "ppm")
## aligner.setParameters(aligner_par)
## aligner.setReference(feature_maps[ref_index])
## 
## for feature_map in feature_maps[:ref_index] + feature_maps[ref_index+1:]:
##     trafo = TransformationDescription() # save the transformed data points
##     aligner.align(feature_map, trafo)
##     trafos[feature_map.getMetaValue("spectra_data")[0].decode()] = trafo
##     transformer = MapAlignmentTransformer()
##     transformer.transformRetentionTimes(feature_map, trafo, True)
## 
## 
## feature_grouper = FeatureGroupingAlgorithmKD()
## 
## consensus_map = ConsensusMap()
## file_descriptions = consensus_map.getColumnHeaders()
## 
## for i, feature_map in enumerate(feature_maps):
##     file_description = file_descriptions.get(i, ColumnHeader())
##     file_description.filename = os.path.basename(
##         feature_map.getMetaValue("spectra_data")[0].decode())
##     file_description.size = feature_map.size()
##     file_descriptions[i] = file_description
## 
## 
## feature_grouper.group(feature_maps, consensus_map)
## consensus_map.setColumnHeaders(file_descriptions)
## consensus_map.setUniqueIds()
## ConsensusXMLFile().store("FeatureMatrix.consensusXML", consensus_map)
## df = consensus_map.get_df()
## df.to_csv('sim3openmsnormal.csv')
## 
## mzML_files = [file3+x for x in os.listdir(file3)]
## 
## feature_maps = []
## for file in mzML_files:
##     # load mzML file into MSExperiment
##     exp = MSExperiment()
##     MzMLFile().load(file, exp) # load each mzML file to an OpenMS file format (MSExperiment)
## 
##     # mass trace detection
##     mass_traces = [] # introduce an empty list where the mass traces will be loaded
##     mtd = MassTraceDetection()
##     mtd_par = mtd.getDefaults() # get the default parameters in order to edit them
##     mtd_par.setValue("mass_error_ppm", 5.0) # high-res instrument, orbitraps
##     mtd_par.setValue("noise_threshold_int", 1.0e02) # data-dependent (usually works for orbitraps)
##     mtd.setParameters(mtd_par) # set the new parameters
##     mtd.run(exp, mass_traces, 0) # run mass trace detection
## 
##     # elution peak detection
##     mass_traces_deconvol = []
##     epd = ElutionPeakDetection()
##     epd_par = epd.getDefaults()
##     epd_par.setValue("width_filtering", "fixed") # The fixed setting filters out mass traces outside the [min_fwhm: 1.0, max_fwhm: 60.0] interval
##     epd.setParameters(epd_par)
##     epd.detectPeaks(mass_traces, mass_traces_deconvol)
## 
##     # feature detection
##     feature_map = FeatureMap() # output features
##     chrom_out = [] # output chromatograms
##     ffm = FeatureFindingMetabo()
##     ffm_par = ffm.getDefaults()
##     #ffm_par.setValue("remove_single_traces", "true") # remove mass traces without satellite isotopic traces
##     #ffm.setParameters(ffm_par)
##     ffm.run(mass_traces_deconvol, feature_map, chrom_out)
##     feature_map.setUniqueIds() # Assigns a new, valid unique id per feature
##     feature_map.setPrimaryMSRunPath([file.encode()]) # Sets the file path to the primary MS run (usually the mzML file)
##     feature_maps.append(feature_map)
## 
## 
## # use as reference for alignment, the file with the largest number of features (works well if you have a pooled QC for example)
## 
## ref_index = feature_maps.index(sorted(feature_maps, key=lambda x: x.size())[-1])
## 
## aligner = MapAlignmentAlgorithmPoseClustering()
## 
## trafos = {}
## 
## # parameter optimization
## aligner_par= aligner.getDefaults()
## aligner_par.setValue("max_num_peaks_considered", -1) # infinite
## aligner_par.setValue("pairfinder:distance_MZ:max_difference", 10.0) # Never pair features with larger m/z distance
## aligner_par.setValue("pairfinder:distance_MZ:unit", "ppm")
## aligner.setParameters(aligner_par)
## aligner.setReference(feature_maps[ref_index])
## 
## for feature_map in feature_maps[:ref_index] + feature_maps[ref_index+1:]:
##     trafo = TransformationDescription() # save the transformed data points
##     aligner.align(feature_map, trafo)
##     trafos[feature_map.getMetaValue("spectra_data")[0].decode()] = trafo
##     transformer = MapAlignmentTransformer()
##     transformer.transformRetentionTimes(feature_map, trafo, True)
## 
## 
## feature_grouper = FeatureGroupingAlgorithmKD()
## 
## consensus_map = ConsensusMap()
## file_descriptions = consensus_map.getColumnHeaders()
## 
## for i, feature_map in enumerate(feature_maps):
##     file_description = file_descriptions.get(i, ColumnHeader())
##     file_description.filename = os.path.basename(
##         feature_map.getMetaValue("spectra_data")[0].decode())
##     file_description.size = feature_map.size()
##     file_descriptions[i] = file_description
## 
## 
## feature_grouper.group(feature_maps, consensus_map)
## consensus_map.setColumnHeaders(file_descriptions)
## consensus_map.setUniqueIds()
## ConsensusXMLFile().store("FeatureMatrix.consensusXML", consensus_map)
## df = consensus_map.get_df()
## df.to_csv('sim3openmstailing.csv')

## ----------------------------------------------------------------------------------------------------
library(mzrtsim)
# the following code will show database in mzrtsim
# MoNA MS1 peaks
data("monams1")
name <- sapply(monams1,function(x) x$name)
length(unique(name))
# MoNA MS1 peaks collected from high resolution mass spectrometry
data("monahrms1")
name <- sapply(monahrms1,function(x) x$name)
length(unique(name))
# HMDB experiment data from GCMS
data("hmdbcms")
name <- sapply(hmdbcms,function(x) x$name)
length(unique(name))
# peak number for different database
pn <- sapply(monams1,function(x) x$np)
mean(as.numeric(pn))
median(as.numeric(pn))
pn <- sapply(hmdbcms,function(x) x$np)
mean(as.numeric(pn))
median(as.numeric(pn))


## ----------------------------------------------------------------------------------------------------
# load detected peaks from xcms, openms, and mzmine
xcms <- read.csv('simxcms.csv')
openms <- read.csv('simopenms.csv')
mzmine <- read.csv('simmzmine.csv')
# load simulated peaks
real <- read.csv('simcsv/case1.csv')
# check overlap
xcmsalign <- enviGCMS::getalign(real$mz,xcms$mz,real$rt,xcms$rt,ppm = 5,deltart = 5)
openmsalign <- enviGCMS::getalign(real$mz,openms$mz,real$rt,openms$RT,ppm = 5,deltart = 5)
mzminealign <- enviGCMS::getalign(real$mz,mzmine$mz,real$rt,mzmine$rt*60,ppm = 5,deltart = 5)
# check unique peaks
xcmsname <- paste(xcmsalign$mz2,xcmsalign$rt2)
openmsname <- paste(openmsalign$mz2,openmsalign$rt2)
mzminename <- paste(mzminealign$mz2,mzminealign$rt2)

length(unique(xcmsname))
# 327/340
length(unique(openmsname))
# 474/564
length(unique(mzminename))
# 523/523

xcmsnamer <- paste(xcmsalign$mz1,xcmsalign$rt1)
openmsnamer <- paste(openmsalign$mz1,openmsalign$rt1)
mzminenamer <- paste(mzminealign$mz1,mzminealign$rt1)

length(unique(xcmsnamer))
# 328/593
length(unique(openmsnamer))
# 449/593
length(unique(mzminenamer))
# 484/593

library(ggvenn)
# display overlap
cvenn <- ggvenn(list(XCMS=real$name[unique(xcmsalign$xid)],OpenMS=real$name[unique(openmsalign$xid)],MZmine4.5=real$name[unique(mzminealign$xid)]))+ggtitle('B')
name <- paste(real$mz,real$rt)
pvenn <- ggvenn(list(XCMS=name[unique(xcmsalign$xid)],OpenMS=name[unique(openmsalign$xid)],MZmine4.5=name[unique(mzminealign$xid)]))+ggtitle('A')

z <- real[real$name %in% unique(real$name[unique(xcmsalign$xid)]),]
library(ggplot2)

data <- data.frame(
  Value = c(real$ins[unique(xcmsalign$xid)], real$ins[-unique(xcmsalign$xid)], real$ins[unique(openmsalign$xid)], real$ins[-unique(openmsalign$xid)], real$ins[unique(mzminealign$xid)], real$ins[-unique(mzminealign$xid)]),
  Peaks = factor(c(rep('TP(XCMS)',length(real$ins[unique(xcmsalign$xid)])),rep('FN(XCMS)',length(real$ins[-unique(xcmsalign$xid)])),rep('TP(OpenMS)',length(real$ins[unique(openmsalign$xid)])),rep('FN(OpenMS)',length(real$ins[-unique(openmsalign$xid)])),rep('TP(MZmine 4.5)',length(real$ins[unique(mzminealign$xid)])),rep('FN(MZmine 4.5)',length(real$ins[-unique(mzminealign$xid)]))))
)

des <- ggplot(data, aes(x = Value, fill = Peaks, color = Peaks)) +
  geom_density(alpha = 0.15) +  
  scale_fill_brewer(palette = "Set1") + 
  scale_color_brewer(palette = "Set1") +
  labs(title = "Overlaid Density Plot", x = "Relative Intensity Distribution", y = "Density") + ggtitle('B') +
  theme_minimal() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

library(patchwork)
p <- pvenn|des
ggsave('figure1.png',p,width = 10,height = 6)



## ----------------------------------------------------------------------------------------------------
# load detected peaks from xcms, openms, and mzmine
xcms <- read.csv('sim3xcmsnormal.csv')
openms <- read.csv('sim3openmsnormal.csv')
mzmine <- read.csv('sim3mzminenormal.csv')
# load simulated peaks
real <- read.csv('simcsv3/normal1.csv')
# check overlap
xcmsalign <- enviGCMS::getalign(real$mz,xcms$mz,real$rt,xcms$rt,ppm = 5,deltart = 5)
openmsalign <- enviGCMS::getalign(real$mz,openms$mz,real$rt,openms$RT,ppm = 5,deltart = 5)
mzminealign <- enviGCMS::getalign(real$mz,mzmine$mz,real$rt,mzmine$rt*60,ppm = 5,deltart = 5)
# check unique peaks
xcmsname <- paste(xcmsalign$mz2,xcmsalign$rt2)
openmsname <- paste(openmsalign$mz2,openmsalign$rt2)
mzminename <- paste(mzminealign$mz2,mzminealign$rt2)

length(unique(xcmsname))
# 367/377
length(unique(openmsname))
# 466/580
length(unique(mzminename))
# 514/516

xcmsnamer <- paste(xcmsalign$mz1,xcmsalign$rt1)
openmsnamer <- paste(openmsalign$mz1,openmsalign$rt1)
mzminenamer <- paste(mzminealign$mz1,mzminealign$rt1)

length(unique(xcmsnamer))
# 367/593
length(unique(openmsnamer))
# 449/593
length(unique(mzminenamer))
# 480/593
ggvenn(list(XCMS=real$name[unique(xcmsalign$xid)],OpenMS=real$name[unique(openmsalign$xid)],MZmine4.5=real$name[unique(mzminealign$xid)]))+ggtitle('B')



## ----------------------------------------------------------------------------------------------------
# load detected peaks from xcms, openms, and mzmine
xcms <- read.csv('sim3xcmsleading.csv')
openms <- read.csv('sim3openmsleading.csv')
mzmine <- read.csv('sim3mzmineleading.csv')
# load simulated peaks
real <- read.csv('simcsv3/leading1.csv')
# check overlap
xcmsalign <- enviGCMS::getalign(real$mz,xcms$mz,real$rt,xcms$rt,ppm = 5,deltart = 5)
openmsalign <- enviGCMS::getalign(real$mz,openms$mz,real$rt,openms$RT,ppm = 5,deltart = 5)
mzminealign <- enviGCMS::getalign(real$mz,mzmine$mz,real$rt,mzmine$rt*60,ppm = 5,deltart = 5)
# check unique peaks
xcmsname <- paste(xcmsalign$mz2,xcmsalign$rt2)
openmsname <- paste(openmsalign$mz2,openmsalign$rt2)
mzminename <- paste(mzminealign$mz2,mzminealign$rt2)

length(unique(xcmsname))
# 396/412
length(unique(openmsname))
# 466/795
length(unique(mzminename))
# 503/503

xcmsnamer <- paste(xcmsalign$mz1,xcmsalign$rt1)
openmsnamer <- paste(openmsalign$mz1,openmsalign$rt1)
mzminenamer <- paste(mzminealign$mz1,mzminealign$rt1)

length(unique(xcmsnamer))
# 397/593
length(unique(openmsnamer))
# 445/593
length(unique(mzminenamer))
# 482/593
ggvenn(list(XCMS=real$name[unique(xcmsalign$xid)],OpenMS=real$name[unique(openmsalign$xid)],MZmine3=real$name[unique(mzminealign$xid)]))+ggtitle('B')



## ----------------------------------------------------------------------------------------------------
# load detected peaks from xcms, openms, and mzmine
xcms <- read.csv('sim3xcmstailing.csv')
openms <- read.csv('sim3openmstailing.csv')
mzmine <- read.csv('sim3mzminetailing.csv')
# load simulated peaks
real <- read.csv('simcsv3/tailing1.csv')
# check overlap
xcmsalign <- enviGCMS::getalign(real$mz,xcms$mz,real$rt,xcms$rt,ppm = 5,deltart = 5)
openmsalign <- enviGCMS::getalign(real$mz,openms$mz,real$rt,openms$RT,ppm = 5,deltart = 5)
mzminealign <- enviGCMS::getalign(real$mz,mzmine$mz,real$rt,mzmine$rt*60,ppm = 5,deltart = 5)
# check unique peaks
xcmsname <- paste(xcmsalign$mz2,xcmsalign$rt2)
openmsname <- paste(openmsalign$mz2,openmsalign$rt2)
mzminename <- paste(mzminealign$mz2,mzminealign$rt2)

length(unique(xcmsname))
# 289/294
length(unique(openmsname))
# 462/567
length(unique(mzminename))
# 511/512

xcmsnamer <- paste(xcmsalign$mz1,xcmsalign$rt1)
openmsnamer <- paste(openmsalign$mz1,openmsalign$rt1)
mzminenamer <- paste(mzminealign$mz1,mzminealign$rt1)

length(unique(xcmsnamer))
# 289/593
length(unique(openmsnamer))
# 443/593
length(unique(mzminenamer))
# 479/593
ggvenn(list(XCMS=real$name[unique(xcmsalign$xid)],OpenMS=real$name[unique(openmsalign$xid)],MZmine3=real$name[unique(mzminealign$xid)]))+ggtitle('B')



## ----------------------------------------------------------------------------------------------------
real <- read.csv('simcsv2/control10.csv')
realsub <- real[real$rt>=rt[21]&real$rt<=rt[30]|real$rt>=rt[51]&real$rt<=rt[70],]
length(unique(realsub$name))
# 28

full <- read.csv('sim2xcms.csv')
cutoff <- read.csv('simxcms.csv')

library(genefilter)
rda <- rowttests(as.matrix(full[,-c(1:3)]),fac=as.factor(c(rep('case',10),rep('control',10))))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05,na.rm = T)
# 916
rt <- seq(10,590,length.out=100)

fullc <- full[which(p.value<0.05),]
sum((fullc$rtmed>=rt[21]&fullc$rtmed<=rt[30])|(fullc$rtmed>=rt[51]&fullc$rtmed<=rt[70]),na.rm = T)
# 781/916
fullchange <- full[full$rtmed>=rt[21]&full$rtmed<=rt[30]|full$rtmed>=rt[51]&full$rtmed<=rt[70],]
sum(full$rtmed>=rt[21]&full$rtmed<=rt[30]|full$rtmed>=rt[51]&full$rtmed<=rt[70])
# 785

align <- enviGCMS::getalign(real$mz,full$mz,real$rt,full$rt,ppm = 5,deltart = 5)
length(unique(real$name[unique(align$xid)]))
# 90
alignc <- enviGCMS::getalign(real$mz,fullc$mzmed,real$rt,fullc$rtmed,ppm = 5,deltart = 5)
length(unique(real$name[unique(alignc$xid)]))
# 43
alignchange <- enviGCMS::getalign(real$mz,fullchange$mzmed,real$rt,fullchange$rtmed,ppm = 5,deltart = 5)
length(unique(paste(fullchange$mz,fullchange$rt)))
# 785
length(unique(real$name[unique(alignchange$xid)]))
# 24
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 24 True positive 19 false positive 4 false negative

real <- read.csv('simcsv/case1.csv')
realsub <- real[real$rt>=rt[21]&real$rt<=rt[30]|real$rt>=rt[51]&real$rt<=rt[70],]
length(unique(realsub$name))
# 28

rda <- rowttests(as.matrix(cutoff[,-c(1:3)]),fac=as.factor(c(rep('case',10),rep('control',10))))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05,na.rm = T)
# 134
cutoffc <- cutoff[which(p.value<0.05),]
sum((cutoffc$rtmed>=rt[21]&cutoffc$rtmed<=rt[30])|(cutoffc$rtmed>=rt[51]&cutoffc$rtmed<=rt[70]),na.rm = T)
# 109
cutoffchange <- cutoff[cutoff$rtmed>=rt[21]&cutoff$rtmed<=rt[30]|cutoff$rtmed>=rt[51]&cutoff$rtmed<=rt[70],]

align <- enviGCMS::getalign(real$mz,cutoff$mzmed,real$rt,cutoff$rtmed,ppm = 5,deltart = 5)
length(unique(real$name[unique(align$xid)]))
# 82
alignc <- enviGCMS::getalign(real$mz,cutoffc$mzmed,real$rt,cutoffc$rtmed,ppm = 5,deltart = 5)
length(unique(real$name[unique(alignc$xid)]))
# 24
alignchange <- enviGCMS::getalign(real$mz,cutoffchange$mz,real$rt,cutoffchange$rt,ppm = 5,deltart = 5)
length(unique(paste(cutoffchange$mz,cutoffchange$rt)))
# 109
length(unique(real$name[unique(alignchange$xid)]))
# 22
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 22 True positive 2 false positive 4 false negative


## ----------------------------------------------------------------------------------------------------
full <- read.csv('sim2openms.csv')
cutoff <- read.csv('simopenms.csv')
library(genefilter)
rda <- rowttests(as.matrix(full[,-c(1:6)]),fac=as.factor(sapply(strsplit(colnames(full[,-c(1:6)]),'[0123456789]'),function(x) x[1])))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05,na.rm = T)
# 2022
rt <- seq(10,590,length.out=100)
fullc <- full[which(p.value<0.05),]
sum((fullc$RT>=rt[21]&fullc$RT<=rt[30])|(fullc$RT>=rt[51]&fullc$RT<=rt[70]))
# 1378
fullchange <- full[full$RT>=rt[21]&full$RT<=rt[30]|full$RT>=rt[51]&full$RT<=rt[70],]
length(unique(paste(fullchange$mz,fullchange$rt)))
# 1387
# check compounds level
real <- read.csv('simcsv2/control10.csv')
align <- enviGCMS::getalign(real$mz,full$mz,real$rt,full$RT,ppm = 5,deltart = 5)
length(unique(real$name[unique(align$xid)]))
# 99
alignc <- enviGCMS::getalign(real$mz,fullc$mz,real$rt,fullc$RT,ppm = 5,deltart = 5)
length(unique(real$name[unique(alignc$xid)]))
# 59 detected changed compounds
alignchange <- enviGCMS::getalign(real$mz,fullchange$mz,real$rt,fullchange$RT,ppm = 5,deltart = 5)
length(unique(real$name[unique(alignchange$xid)]))
# 29 changed compounds
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 28 True positive 31 False positive 

rda <- rowttests(as.matrix(cutoff[,-c(1:6)]),fac=as.factor(sapply(strsplit(colnames(cutoff[,-c(1:6)]),'[0123456789]'),function(x) x[1])))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05,na.rm = T)
# 201
cutoffc <- cutoff[which(p.value<0.05),]
sum((cutoffc$RT>=rt[21]&cutoffc$RT<=rt[30])|(cutoffc$RT>=rt[51]&cutoffc$RT<=rt[70]))
# 162
cutoffchange <- cutoff[cutoff$RT>=rt[21]&cutoff$RT<=rt[30]|cutoff$RT>=rt[51]&cutoff$RT<=rt[70],]
length(unique(paste(cutoffchange$mz,cutoffchange$rt)))
# 169
real <- read.csv('simcsv/control10.csv')
align <- enviGCMS::getalign(real$mz,cutoff$mz,real$rt,cutoff$RT,ppm = 5,deltart = 5)
length(unique(real$name[unique(align$xid)]))
# 99 
alignc <- enviGCMS::getalign(real$mz,cutoffc$mz,real$rt,cutoffc$RT,ppm = 5,deltart = 5)
length(unique(real$name[unique(alignc$xid)]))
# 28 detected changed compounds
alignchange <- enviGCMS::getalign(real$mz,cutoffchange$mz,real$rt,cutoffchange$RT,ppm = 5,deltart = 5)
length(unique(real$name[unique(alignchange$xid)]))
# 28 changed compounds
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 28 True positive 1 False positive


## ----------------------------------------------------------------------------------------------------
full <- read.csv('sim2mzmine.csv')
cutoff <- read.csv('simmzmine.csv')

fulldata <- full[,grepl('datafile(.*?)height',colnames(full))]
cutoffdata <- cutoff[,grepl('datafile(.*?)height',colnames(full))]

library(genefilter)
rda <- rowttests(as.matrix(fulldata),fac=as.factor(sapply(strsplit(colnames(fulldata),'\\.|[0123456789]'),function(x) x[2])))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05,na.rm = T)
# 1484
rt <- seq(10,590,length.out=100)

fullc <- full[which(p.value<0.05),]
sum((fullc$rt*60>=rt[21]&fullc$rt*60<=rt[30])|(fullc$rt*60>=rt[51]&fullc$rt*60<=rt[70]))
# 1404
fullchange <- full[full$rt*60>=rt[21]&full$rt*60<=rt[30]|full$rt*60>=rt[51]&full$rt*60<=rt[70],]
length(unique(paste(fullchange$mz,fullchange$rt)))
# 1555

# check compounds level
real <- read.csv('simcsv2/control10.csv')
align <- enviGCMS::getalign(real$mz,full$mz,real$rt,full$rt*60,ppm = 5,deltart = 5)
length(unique(real$name[unique(align$xid)]))
# 99
alignc <- enviGCMS::getalign(real$mz,fullc$mz,real$rt,fullc$rt*60,ppm = 5,deltart = 5)
length(unique(real$name[unique(alignc$xid)]))
# 42 detected changed compounds
alignchange <- enviGCMS::getalign(real$mz,fullchange$mz,real$rt,fullchange$rt*60,ppm = 5,deltart = 5)
length(unique(real$name[unique(alignchange$xid)]))
# 28 changed compounds
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 28 True positive 14 False positive

rda <- rowttests(as.matrix(cutoffdata),fac=as.factor(sapply(strsplit(colnames(cutoffdata),'\\.|[0123456789]'),function(x) x[2])))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05,na.rm = T)
# 208
rt <- seq(10,590,length.out=100)
cutoffc <- cutoff[which(p.value<0.05),]
sum((cutoffc$rt*60>=rt[21]&cutoffc$rt*60<=rt[30])|(cutoffc$rt*60>=rt[51]&cutoffc$rt*60<=rt[70]))
# 202
cutoffchange <- cutoff[cutoff$rt*60>=rt[21]&cutoff$rt*60<=rt[30]|cutoff$rt*60>=rt[51]&cutoff$rt*60<=rt[70],]
length(unique(paste(cutoffchange$mz,fullchange$rt)))
# 1161

# check compounds level
real <- read.csv('simcsv2/control10.csv')
align <- enviGCMS::getalign(real$mz,cutoff$mz,real$rt,cutoff$rt*60,ppm = 5,deltart = 5)
length(unique(real$name[unique(align$xid)]))
# 99
alignc <- enviGCMS::getalign(real$mz,cutoffc$mz,real$rt,cutoffc$rt*60,ppm = 5,deltart = 5)
length(unique(real$name[unique(alignc$xid)]))
# 29 detected changed compounds
alignchange <- enviGCMS::getalign(real$mz,cutoffchange$mz,real$rt,cutoffchange$rt*60,ppm = 5,deltart = 5)
length(unique(real$name[unique(alignchange$xid)]))
# 28 changed compounds
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 28 True positive 1 false positive


## ----------------------------------------------------------------------------------------------------
knitr::knit_exit()


## ----------------------------------------------------------------------------------------------------
library(mzrtsim)
# load the high resolution MS1 database from MoNA
data(monahrms1)
# Set peak height
ph1 <- c(15,15,15)
rt <- c(120,200,250)
compound <- c(1:3)
rf <- c(122,222,345)
simmzml(name='toc',db=monahrms1,pheight = ph1,compound=compound,rtime = rt, rf=rf,unique = T,matrix = T,tailingfactor = c(0.5,1,1.5))
mzrange=c(500,600)
rtrange=c(100,300)
mzml = 'sim/case/case1.mzML'
dt <- Spectra::Spectra(mzml)
rt <- Spectra::rtime(dt)
ins <- Spectra::intensity(dt)
mz <- Spectra::mz(dt)
mzv <- unlist(mz)
rtimev <- rep(rt, times = sapply(mz, length))
intensityv <- log(unlist(ins)+1)
idx <- mzv>200&mzv<400
mzv <- mzv[idx]
rtimev <- rtimev[idx]
intensityv <- intensityv[idx]
norm <- (intensityv - min(intensityv)) / (max(intensityv) - min(intensityv))
plot(
                        rtimev,
                        mzv,
                        pch = 15,
                        cex = 0.1,
                        col = grDevices::gray(1 - norm),
                        xlab = 'retention time(s)',
                        ylab = 'm/z'
                )                


## ----------------------------------------------------------------------------------------------------
realm <- real[-xcmsalign$xid,]
realn <- real[xcmsalign$xid,]
xcmsalign2 <- enviGCMS::getalign(real$mz,xcms$mz,real$rt,xcms$rt)
xcmsm <- enviGCMS::getfilter(xcms,unique(xcmsalign2$xid))
xcmsm2 <- enviGCMS::getfilter(xcms,-unique(xcmsalign2$xid))
sin <- apply(xcmsm$data,1,mean)
sin2 <- apply(xcmsm2$data,1,mean)
set.seed(1)
compound <- sample(c(1:1114),100)
library(mzrtsim)
data("monahrms1")
uniquecpidx <- sapply(monahrms1, function(x) x$name)
monahrms1 <- monahrms1[!duplicated(uniquecpidx)]
x <- monahrms1[compound]
sl <- sapply(x,function(x) x$spectra)
mz <- lapply(x,function(x) x$spectra$mz)
ins <- lapply(x, function(x) x$spectra$ins)
mzl <- do.call(c,mz)
insl <- do.call(c,ins)
rt <- seq(10,590,length.out=100)
length <- lapply(x,function(x) length(x$spectra$mz))
rtt <- rep(rt,sapply(mz,length))
raw <- cbind.data.frame(rt=rtt,mz=mzl,ins=insl)
rawd <- raw[raw$ins>5&raw$mz>100&raw$mz<1000,]

rawx<- enviGCMS::getalign(rawd$mz,realm$mz,rawd$rt,realm$rt)
rawy<- enviGCMS::getalign(rawd$mz,realn$mz,rawd$rt,realn$rt)
z <- rawd[unique(rawx$xid),]
zz <- rawd[unique(rawy$xid),]


## ----------------------------------------------------------------------------------------------------
full <- read.csv('sim2xcms.csv')
cutoff <- read.csv('simxcms.csv')
# impute NA just in case
full[is.na(full)] <- min(full[,-c(1:3)],na.rm = T)
cutoff[is.na(cutoff)] <- min(cutoff[,-c(1:3)],na.rm = T)

library(genefilter)
rda <- rowttests(as.matrix(full[,-c(1:3)]),fac=as.factor(c(rep('case',10),rep('control',10))))
p.value <- p.adjust(rda$p.value,'bonferroni')
sum(p.value<0.05,na.rm = T)
# 917
rt <- seq(10,590,length.out=100)

fullc <- full[which(p.value<0.05),]
sum((fullc$rtmed>=rt[21]&fullc$rtmed<=rt[30])|(fullc$rtmed>=rt[51]&fullc$rtmed<=rt[70]))
# 782
fullchange <- full[full$rtmed>=rt[21]&full$rtmed<=rt[30]|full$rtmed>=rt[51]&full$rtmed<=rt[70],]

real <- read.csv('simcsv2/control10.csv')
align <- enviGCMS::getalign(real$mz,full$mz,real$rt,full$rt,ppm = 5,deltart = 5)
length(unique(real$name[unique(align$xid)]))
# 90
alignc <- enviGCMS::getalign(real$mz,fullc$mzmed,real$rt,fullc$rtmed,ppm = 5,deltart = 5)
length(unique(real$name[unique(alignc$xid)]))
# 43
alignchange <- enviGCMS::getalign(real$mz,fullchange$mzmed,real$rt,fullchange$rtmed,ppm = 5,deltart = 5)
length(unique(paste(fullchange$mz,fullchange$rt)))
# 785
length(unique(real$name[unique(alignchange$xid)]))
# 24
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 24 True positive 19 false positive  

rda <- rowttests(as.matrix(cutoff[,-c(1:3)]),fac=as.factor(c(rep('case',10),rep('control',10))))
p.value <- p.adjust(rda$p.value,'bonferroni')
sum(p.value<0.05,na.rm = T)
# 134
cutoffc <- cutoff[which(p.value<0.05),]
sum((cutoffc$rtmed>=rt[21]&cutoffc$rtmed<=rt[30])|(cutoffc$rtmed>=rt[51]&cutoffc$rtmed<=rt[70]))
# 109
cutoffchange <- cutoff[cutoff$rtmed>=rt[21]&cutoff$rtmed<=rt[30]|cutoff$rtmed>=rt[51]&cutoff$rtmed<=rt[70],]

real <- read.csv('simcsv/case1.csv')
align <- enviGCMS::getalign(real$mz,cutoff$mzmed,real$rt,cutoff$rtmed,ppm = 5,deltart = 5)
length(unique(real$name[unique(align$xid)]))
# 82
alignc <- enviGCMS::getalign(real$mz,cutoffc$mzmed,real$rt,cutoffc$rtmed,ppm = 5,deltart = 5)
length(unique(real$name[unique(alignc$xid)]))
# 24
alignchange <- enviGCMS::getalign(real$mz,cutoffchange$mz,real$rt,cutoffchange$rt,ppm = 5,deltart = 5)
length(unique(paste(cutoffchange$mz,cutoffchange$rt)))
# 109
length(unique(real$name[unique(alignchange$xid)]))
# 22
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 22 True positive 2 false positive
