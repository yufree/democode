## ----setup, include=FALSE-----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## -----------------------------------------------------------------------------------------------------
library(mzrtsim)
# load the high resolution MS1 database from MoNA
data(monahrms1)
# Set peak width
pw1 <- c(rep(5,30),rep(10,40),rep(15,30))
pw2 <- c(rep(5,20),rep(10,30),rep(15,50))
# Set retention time
rt <- seq(10,590,length.out=100)
# rt[c(21:30,51:70)]
# display peaks profile
plot(pw1~rt,cex=0.5,col='blue')
points(pw2~rt,col='red')
# for reproducible purpose
set.seed(1)
# select compounds
compound <- sample(c(1:1114),100)
set.seed(2)
# sample s/n for each compounds
SNR <- sample(c(100:10000),100)

# 30% changed
for(i in c(1:10)){
  # with unique, only one spectra will be used for simulated compound and no duplicated spectra for the same compound as database will contain multiple spectra for the same compound
  simmzml(name=paste0('sim/case/case',i),db=monahrms1,pwidth = pw1,compound=compound,rtime = rt, SNR=SNR,unique = T)
}

for(i in c(1:10)){
  # set different peak width
  simmzml(name=paste0('sim/control/control',i),db=monahrms1,pwidth = pw2,compound=compound,rtime = rt, SNR=SNR,unique = T)
}


## -----------------------------------------------------------------------------------------------------
library(mzrtsim)
# the following part is the same with previous simulation
data(monahrms1)
pw1 <- c(rep(5,30),rep(10,40),rep(15,30))
pw2 <- c(rep(5,20),rep(10,30),rep(15,50))
rt <- seq(10,590,length.out=100)
rt[c(21:30,51:70)]
plot(pw1~rt,cex=0.5,col='blue')
points(pw2~rt,col='red')
set.seed(1)
compound <- sample(c(1:1114),100)
set.seed(2)
SNR <- sample(c(100:10000),100)

# 30% changed
for(i in c(1:10)){
  # set intenisty cutoff as 0
  simmzml(name=paste0('sim4/case/case',i),db=monahrms1,pwidth = pw1,compound=compound,rtime = rt, SNR=SNR,unique = T,inscutoff = 0)
}

for(i in c(1:10)){
  # set intenisty cutoff as 0
  simmzml(name=paste0('sim4/control/control',i),db=monahrms1,pwidth = pw2,compound=compound,rtime = rt, SNR=SNR,unique = T,inscutoff = 0)
}

files <- list.files('sim4',full.names = T,recursive = T,pattern = '*.csv')
name <- basename(files)
file.rename(files,paste0('simcsv4/',name))


## -----------------------------------------------------------------------------------------------------
library(mzrtsim)
# the following part is the same with previous simulation
data(monahrms1)
pw1 <- c(rep(5,30),rep(10,40),rep(15,30))
pw2 <- c(rep(5,30),rep(10,40),rep(15,30))
pw3 <- c(rep(5,30),rep(10,40),rep(15,30))
rt <- seq(10,590,length.out=100)
plot(pw1~rt,cex=0.5,col='blue')
points(pw2~rt,col='red')
set.seed(1)
compound <- sample(c(1:1114),100)
set.seed(2)
SNR <- sample(c(100:10000),100)

for(i in c(1:10)){
  # change tailing factor to 1.5 to simulate tailing peaks, add matrix effect
  simmzml(name=paste0('sim5/tailing/tailing',i),db=monahrms1,pwidth = pw1,compound=compound,rtime = rt, SNR=SNR,unique = T,tailingfactor = 1.5,matrix = T)
}

for(i in c(1:10)){
  # change tailing factor to 1 to simulate normal peaks, add matrix effect
  simmzml(name=paste0('sim5/normal/normal',i),db=monahrms1,pwidth = pw2,compound=compound,rtime = rt, SNR=SNR,unique = T,tailingfactor = 1,matrix = T)
}

for(i in c(1:10)){
  # change tailing factor to 0.8 to simulate leading peaks, add matrix effect
  simmzml(name=paste0('sim5/leading/leading',i),db=monahrms1,pwidth = pw2,compound=compound,rtime = rt, SNR=SNR,unique = T,tailingfactor = 0.8,matrix = T)
}

files <- list.files('sim5',full.names = T,recursive = T,pattern = '*.csv')
name <- basename(files)
file.rename(files,paste0('simcsv5/',name))


## -----------------------------------------------------------------------------------------------------
library(xcms)
library(MsExperiment)
# add this for mac os
# setMSnbaseFastLoad(opt = F)
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
# save peaks profile with name
compound <- read.csv('simcsv/case1.csv')
# 593 peaks
# simulate mass range [100,1000] for MS1 scan
compoundsub <- compound[compound$mz>100&compound$mz<1000,]
# 533 peaks
# align simulated peaks with detected peaks
align2 <- enviGCMS::getalign(compoundsub$mz,re$mzmed,compoundsub$rt,re$rtmed)
length(unique(compound$name[align2$xid]))
# 81
# 452 peaks
# 446 peaks
library(SummarizedExperiment)
res <- quantify(xdata, value = "into", method = "sum")
data <- assay(res)
library(genefilter)
# check changed peaks
rda <- rowttests(as.matrix(data),fac=as.factor(pd$sample_group))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05)
df <- cbind.data.frame(re,data)
write.csv(df,'simxcms.csv')


## -----------------------------------------------------------------------------------------------------
library(xcms)
library(MsExperiment)
# add this for mac os
# setMSnbaseFastLoad(opt = F)
# peak picking for simulated peaks
path <- 'sim4'
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
# save peaks profile with name
compound <- read.csv('simcsv4/case1.csv')
# 593 peaks
# simulate mass range [100,1000] for MS1 scan
compoundsub <- compound[compound$mz>100&compound$mz<1000,]
# 533 peaks
# align simulated peaks with detected peaks
align2 <- enviGCMS::getalign(compoundsub$mz,re$mzmed,compoundsub$rt,re$rtmed)
length(unique(compound$name[align2$xid]))
# 70
# 2844
# 2264 peaks
library(SummarizedExperiment)
res <- quantify(xdata, value = "into", method = "sum")
data <- assay(res)
library(genefilter)
# check changed peaks
rda <- rowttests(as.matrix(data),fac=as.factor(pd$sample_group))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05)
# 495
df <- cbind.data.frame(re,data)
write.csv(df,'sim4xcms.csv')


## ----leading------------------------------------------------------------------------------------------
library(xcms)
library(MsExperiment)
# add this for mac os
# setMSnbaseFastLoad(opt = F)
# peak picking for simulated peaks
path <- 'sim5/leading/'
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
# save peaks profile with name
compound <- read.csv('simcsv5/leading1.csv')
# 593 peaks
# simulate mass range [100,1000] for MS1 scan
compoundsub <- compound[compound$mz>100&compound$mz<1000,]
# 533 peaks
# align simulated peaks with detected peaks
align2 <- enviGCMS::getalign(compoundsub$mz,re$mzmed,compoundsub$rt,re$rtmed)
length(unique(compound$name[align2$xid]))
# 83
# 481 peaks
# 484 peaks
library(SummarizedExperiment)
res <- quantify(xdata, value = "into", method = "sum")
data <- assay(res)
leading <- cbind.data.frame(re,data)
write.csv(leading,'sim5xcmsleading.csv')


## ----normal-------------------------------------------------------------------------------------------
library(xcms)
library(MsExperiment)
# add this for mac os
# setMSnbaseFastLoad(opt = F)
# peak picking for simulated peaks
path <- 'sim5/normal'
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
# save peaks profile with name
compound <- read.csv('simcsv5/normal1.csv')
# 593 peaks
# simulate mass range [100,1000] for MS1 scan
compoundsub <- compound[compound$mz>100&compound$mz<1000,]
# 533 peaks
# align simulated peaks with detected peaks
align2 <- enviGCMS::getalign(compoundsub$mz,re$mzmed,compoundsub$rt,re$rtmed)
length(unique(compoundsub$name[align2$xid]))
# 85
# 453 peaks
# 458 peaks
library(SummarizedExperiment)
res <- quantify(xdata, value = "into", method = "sum")
data <- assay(res)
normal <- cbind.data.frame(re,data)
write.csv(normal,'sim5xcmsnormal.csv')


## ----tailing------------------------------------------------------------------------------------------
library(xcms)
library(MsExperiment)
# add this for mac os
# setMSnbaseFastLoad(opt = F)
# peak picking for simulated peaks
path <- 'sim5/tailing'
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
# save peaks profile with name
compound <- read.csv('simcsv5/tailing1.csv')
# 593 peaks
# simulate mass range [100,1000] for MS1 scan
compoundsub <- compound[compound$mz>100&compound$mz<1000,]
# 533 peaks
# align simulated peaks with detected peaks
align2 <- enviGCMS::getalign(compoundsub$mz,re$mzmed,compoundsub$rt,re$rtmed)
length(unique(compound$name[align2$xid]))
# 71
# 357 peaks
# 355 peaks
library(SummarizedExperiment)
res <- quantify(xdata, value = "into", method = "sum")
data <- assay(res)
tailing <- cbind.data.frame(re,data)
write.csv(tailing,'sim5xcmstailing.csv')


## -----------------------------------------------------------------------------------------------------
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
## file1 = "sim4/case/"
## file2 = "sim4/control/"
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
## df.to_csv('sim4openms.csv')

## import os
## import shutil
## import requests
## import pandas as pd
## from pyopenms import *
## import os
## import numpy as np
## file1 = "sim5/leading/"
## file2 = "sim5/normal/"
## file3 = "sim5/tailing/"
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
## df.to_csv('sim5openmsleading.csv')
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
## df.to_csv('sim5openmsnormal.csv')
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
## df.to_csv('sim5openmstailing.csv')

## -----------------------------------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------------------------------
# load detected peaks from xcms, openms, and mzmine
xcms <- read.csv('simxcms.csv')
openms <- read.csv('simopenms.csv')
mzmine <- read.csv('simmzmine.csv')
# load simulated peaks
real <- read.csv('simcsv/case1.csv')
# check overlap
xcmsalign <- enviGCMS::getalign(real$mz,xcms$mz,real$rt,xcms$rt)
openmsalign <- enviGCMS::getalign(real$mz,openms$mz,real$rt,openms$RT)
mzminealign <- enviGCMS::getalign(real$mz,mzmine$mz,real$rt,mzmine$rt*60)
# check unique peaks
xcmsname <- paste(xcmsalign$mz2,xcmsalign$rt2)
openmsname <- paste(openmsalign$mz2,openmsalign$rt2)
mzminename <- paste(mzminealign$mz2,mzminealign$rt2)

length(unique(xcmsname))
# 443/446
length(unique(openmsname))
# 657/1025
length(unique(mzminename))
# 1199/28222

xcmsnamer <- paste(xcmsalign$mz1,xcmsalign$rt1)
openmsnamer <- paste(openmsalign$mz1,openmsalign$rt1)
mzminenamer <- paste(mzminealign$mz1,mzminealign$rt1)

length(unique(xcmsnamer))
# 446/593
length(unique(openmsnamer))
# 487/593
length(unique(mzminenamer))
# 504/593

library(ggvenn)
# display overlap
cvenn <- ggvenn(list(XCMS=real$name[unique(xcmsalign$xid)],OpenMS=real$name[unique(openmsalign$xid)],MZmine3=real$name[unique(mzminealign$xid)]))+ggtitle('B')
name <- paste(real$mz,real$rt)
pvenn <- ggvenn(list(XCMS=name[unique(xcmsalign$xid)],OpenMS=name[unique(openmsalign$xid)],MZmine3=name[unique(mzminealign$xid)]))+ggtitle('A')

library(ggplot2)

data <- data.frame(
  Value = c(real$ins[unique(xcmsalign$xid)], real$ins[-unique(xcmsalign$xid)], real$ins[unique(openmsalign$xid)], real$ins[-unique(openmsalign$xid)], real$ins[unique(mzminealign$xid)], real$ins[-unique(mzminealign$xid)]),
  Peaks = factor(c(rep('TP(XCMS)',length(real$ins[unique(xcmsalign$xid)])),rep('FN(XCMS)',length(real$ins[-unique(xcmsalign$xid)])),rep('TP(OpenMS)',length(real$ins[unique(openmsalign$xid)])),rep('FN(OpenMS)',length(real$ins[-unique(openmsalign$xid)])),rep('TP(MZmine 3)',length(real$ins[unique(mzminealign$xid)])),rep('FN(MZmine 3)',length(real$ins[-unique(mzminealign$xid)]))))
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



## -----------------------------------------------------------------------------------------------------
# load detected peaks from xcms, openms, and mzmine
xcms <- read.csv('sim5xcmsnormal.csv')
openms <- read.csv('sim5openmsnormal.csv')
mzmine <- read.csv('sim5mzminenormal.csv')
# load simulated peaks
real <- read.csv('simcsv5/normal1.csv')
# check overlap
xcmsalign <- enviGCMS::getalign(real$mz,xcms$mz,real$rt,xcms$rt)
openmsalign <- enviGCMS::getalign(real$mz,openms$mz,real$rt,openms$RT)
mzminealign <- enviGCMS::getalign(real$mz,mzmine$mz,real$rt,mzmine$rt*60)
# check unique peaks
xcmsname <- paste(xcmsalign$mz2,xcmsalign$rt2)
openmsname <- paste(openmsalign$mz2,openmsalign$rt2)
mzminename <- paste(mzminealign$mz2,mzminealign$rt2)

length(unique(xcmsname))
# 448/453
length(unique(openmsname))
# 646/971
length(unique(mzminename))
# 999/15864

xcmsnamer <- paste(xcmsalign$mz1,xcmsalign$rt1)
openmsnamer <- paste(openmsalign$mz1,openmsalign$rt1)
mzminenamer <- paste(mzminealign$mz1,mzminealign$rt1)

length(unique(xcmsnamer))
# 451/593
length(unique(openmsnamer))
# 475/593
length(unique(mzminenamer))
# 504/593
ggvenn(list(XCMS=real$name[unique(xcmsalign$xid)],OpenMS=real$name[unique(openmsalign$xid)],MZmine3=real$name[unique(mzminealign$xid)]))+ggtitle('B')



## -----------------------------------------------------------------------------------------------------
# load detected peaks from xcms, openms, and mzmine
xcms <- read.csv('sim5xcmsleading.csv')
openms <- read.csv('sim5openmsleading.csv')
mzmine <- read.csv('sim5mzmineleading.csv')
# load simulated peaks
real <- read.csv('simcsv5/leading1.csv')
# check overlap
xcmsalign <- enviGCMS::getalign(real$mz,xcms$mz,real$rt,xcms$rt)
openmsalign <- enviGCMS::getalign(real$mz,openms$mz,real$rt,openms$RT)
mzminealign <- enviGCMS::getalign(real$mz,mzmine$mz,real$rt,mzmine$rt*60)
# check unique peaks
xcmsname <- paste(xcmsalign$mz2,xcmsalign$rt2)
openmsname <- paste(openmsalign$mz2,openmsalign$rt2)
mzminename <- paste(mzminealign$mz2,mzminealign$rt2)

length(unique(xcmsname))
# 474/481
length(unique(openmsname))
# 641/1181
length(unique(mzminename))
# 982/15855

xcmsnamer <- paste(xcmsalign$mz1,xcmsalign$rt1)
openmsnamer <- paste(openmsalign$mz1,openmsalign$rt1)
mzminenamer <- paste(mzminealign$mz1,mzminealign$rt1)

length(unique(xcmsnamer))
# 481/593
length(unique(openmsnamer))
# 483/593
length(unique(mzminenamer))
# 505/593
ggvenn(list(XCMS=real$name[unique(xcmsalign$xid)],OpenMS=real$name[unique(openmsalign$xid)],MZmine3=real$name[unique(mzminealign$xid)]))+ggtitle('B')



## -----------------------------------------------------------------------------------------------------
# load detected peaks from xcms, openms, and mzmine
xcms <- read.csv('sim5xcmstailing.csv')
openms <- read.csv('sim5openmstailing.csv')
mzmine <- read.csv('sim5mzminetailing.csv')
# load simulated peaks
real <- read.csv('simcsv5/tailing1.csv')
# check overlap
xcmsalign <- enviGCMS::getalign(real$mz,xcms$mz,real$rt,xcms$rt)
openmsalign <- enviGCMS::getalign(real$mz,openms$mz,real$rt,openms$RT)
mzminealign <- enviGCMS::getalign(real$mz,mzmine$mz,real$rt,mzmine$rt*60)
# check unique peaks
xcmsname <- paste(xcmsalign$mz2,xcmsalign$rt2)
openmsname <- paste(openmsalign$mz2,openmsalign$rt2)
mzminename <- paste(mzminealign$mz2,mzminealign$rt2)

length(unique(xcmsname))
# 350/355
length(unique(openmsname))
# 577/850
length(unique(mzminename))
# 1021/15918

xcmsnamer <- paste(xcmsalign$mz1,xcmsalign$rt1)
openmsnamer <- paste(openmsalign$mz1,openmsalign$rt1)
mzminenamer <- paste(mzminealign$mz1,mzminealign$rt1)

length(unique(xcmsnamer))
# 351/593
length(unique(openmsnamer))
# 471/593
length(unique(mzminenamer))
# 504/593
ggvenn(list(XCMS=real$name[unique(xcmsalign$xid)],OpenMS=real$name[unique(openmsalign$xid)],MZmine3=real$name[unique(mzminealign$xid)]))+ggtitle('B')



## -----------------------------------------------------------------------------------------------------
full <- read.csv('sim4xcms.csv')
cutoff <- read.csv('simxcms.csv')
# impute NA just in case
full[is.na(full)] <- min(full[,-c(1:3)],na.rm = T)
cutoff[is.na(cutoff)] <- min(cutoff[,-c(1:3)],na.rm = T)

library(genefilter)
rda <- rowttests(as.matrix(full[,-c(1:3)]),fac=as.factor(c(rep('case',10),rep('control',10))))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05,na.rm = T)
# 495
rt <- seq(10,590,length.out=100)

fullc <- full[p.value<0.05,]
sum((fullc$rtmed>=rt[21]&fullc$rtmed<=rt[30])|(fullc$rtmed>=rt[51]&fullc$rtmed<=rt[70]))
# 458
fullchange <- full[full$rtmed>=rt[21]&full$rtmed<=rt[30]|full$rtmed>=rt[51]&full$rtmed<=rt[70],]

real <- read.csv('simcsv4/control10.csv')
align <- enviGCMS::getalign(real$mz,full$mz,real$rt,full$rt)
length(unique(real$name[unique(align$xid)]))
# 93
alignc <- enviGCMS::getalign(real$mz,fullc$mzmed,real$rt,fullc$rtmed)
length(unique(real$name[unique(alignc$xid)]))
# 31
alignchange <- enviGCMS::getalign(real$mz,fullchange$mzmed,real$rt,fullchange$rtmed)
length(unique(paste(fullchange$mz,fullchange$rt)))
# 797
length(unique(real$name[unique(alignchange$xid)]))
# 26
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 25 True positive 6 false positive 1 false negative  

rda <- rowttests(as.matrix(cutoff[,-c(1:3)]),fac=as.factor(c(rep('case',10),rep('control',10))))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05)
# 156
cutoffc <- cutoff[p.value<0.05,]
sum((cutoffc$rtmed>=rt[21]&cutoffc$rtmed<=rt[30])|(cutoffc$rtmed>=rt[51]&cutoffc$rtmed<=rt[70]))
# 123
cutoffchange <- cutoff[cutoff$rtmed>=rt[21]&cutoff$rtmed<=rt[30]|cutoff$rtmed>=rt[51]&cutoff$rtmed<=rt[70],]

real <- read.csv('simcsv/case1.csv')
align <- enviGCMS::getalign(real$mz,cutoff$mzmed,real$rt,cutoff$rtmed)
length(unique(real$name[unique(align$xid)]))
# 89
alignc <- enviGCMS::getalign(real$mz,cutoffc$mzmed,real$rt,cutoffc$rtmed)
length(unique(real$name[unique(alignc$xid)]))
# 32
alignchange <- enviGCMS::getalign(real$mz,cutoffchange$mz,real$rt,cutoffchange$rt)
length(unique(paste(cutoffchange$mz,cutoffchange$rt)))
# 145
length(unique(real$name[unique(alignchange$xid)]))
# 27
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 25 True positive 7 false positive 2 false negative


## -----------------------------------------------------------------------------------------------------
full <- read.csv('sim4xcms.csv')
cutoff <- read.csv('simxcms.csv')
# impute NA just in case
full[is.na(full)] <- min(full[,-c(1:3)],na.rm = T)
cutoff[is.na(cutoff)] <- min(cutoff[,-c(1:3)],na.rm = T)

library(genefilter)
rda <- rowttests(as.matrix(full[,-c(1:3)]),fac=as.factor(c(rep('case',10),rep('control',10))))
p.value <- p.adjust(rda$p.value,'bonferroni')
sum(p.value<0.05,na.rm = T)
# 261
rt <- seq(10,590,length.out=100)

fullc <- full[p.value<0.05,]
sum((fullc$rtmed>=rt[21]&fullc$rtmed<=rt[30])|(fullc$rtmed>=rt[51]&fullc$rtmed<=rt[70]))
# 244
fullchange <- full[full$rtmed>=rt[21]&full$rtmed<=rt[30]|full$rtmed>=rt[51]&full$rtmed<=rt[70],]

real <- read.csv('simcsv4/control10.csv')
align <- enviGCMS::getalign(real$mz,full$mz,real$rt,full$rt)
length(unique(real$name[unique(align$xid)]))
# 93
alignc <- enviGCMS::getalign(real$mz,fullc$mzmed,real$rt,fullc$rtmed)
length(unique(real$name[unique(alignc$xid)]))
# 20
alignchange <- enviGCMS::getalign(real$mz,fullchange$mzmed,real$rt,fullchange$rtmed)
length(unique(paste(fullchange$mz,fullchange$rt)))
# 797
length(unique(real$name[unique(alignchange$xid)]))
# 26
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 19 True positive 7 false positive 1 false negative  

rda <- rowttests(as.matrix(cutoff[,-c(1:3)]),fac=as.factor(c(rep('case',10),rep('control',10))))
p.value <- p.adjust(rda$p.value,'bonferroni')
sum(p.value<0.05)
# 111
cutoffc <- cutoff[p.value<0.05,]
sum((cutoffc$rtmed>=rt[21]&cutoffc$rtmed<=rt[30])|(cutoffc$rtmed>=rt[51]&cutoffc$rtmed<=rt[70]))
# 88
cutoffchange <- cutoff[cutoff$rtmed>=rt[21]&cutoff$rtmed<=rt[30]|cutoff$rtmed>=rt[51]&cutoff$rtmed<=rt[70],]

real <- read.csv('simcsv/case1.csv')
align <- enviGCMS::getalign(real$mz,cutoff$mzmed,real$rt,cutoff$rtmed)
length(unique(real$name[unique(align$xid)]))
# 89
alignc <- enviGCMS::getalign(real$mz,cutoffc$mzmed,real$rt,cutoffc$rtmed)
length(unique(real$name[unique(alignc$xid)]))
# 20
alignchange <- enviGCMS::getalign(real$mz,cutoffchange$mz,real$rt,cutoffchange$rt)
length(unique(paste(cutoffchange$mz,cutoffchange$rt)))
# 145
length(unique(real$name[unique(alignchange$xid)]))
# 27
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 18 True positive 1 false positive 9 false negative


## -----------------------------------------------------------------------------------------------------
full <- read.csv('sim4openms.csv')
cutoff <- read.csv('simopenms.csv')
# impute NA just in case
full[is.na(full)] <- min(full[,-c(1:6)],na.rm = T)
cutoff[is.na(cutoff)] <- min(cutoff[,-c(1:6)],na.rm = T)
library(genefilter)
rda <- rowttests(as.matrix(full[,-c(1:6)]),fac=as.factor(sapply(strsplit(colnames(full[,-c(1:6)]),'[0123456789]'),function(x) x[1])))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05,na.rm = T)
# 1111
rt <- seq(10,590,length.out=100)
fullc <- full[p.value<0.05,]
sum((fullc$RT>=rt[21]&fullc$RT<=rt[30])|(fullc$RT>=rt[51]&fullc$RT<=rt[70]))
# 1060
fullchange <- full[full$RT>=rt[21]&full$RT<=rt[30]|full$RT>=rt[51]&full$RT<=rt[70],]
length(unique(paste(fullchange$mz,fullchange$rt)))
# 2129
# check compounds level
real <- read.csv('simcsv4/control10.csv')
align <- enviGCMS::getalign(real$mz,full$mz,real$rt,full$RT)
length(unique(real$name[unique(align$xid)]))
# 99
alignc <- enviGCMS::getalign(real$mz,fullc$mz,real$rt,fullc$RT)
length(unique(real$name[unique(alignc$xid)]))
# 35 detected changed compounds
alignchange <- enviGCMS::getalign(real$mz,fullchange$mz,real$rt,fullchange$RT)
length(unique(real$name[unique(alignchange$xid)]))
# 31 changed compounds
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 28 True positive 7 False positive 3 False negative
unique(real$name[unique(alignchange$xid)])

rda <- rowttests(as.matrix(cutoff[,-c(1:6)]),fac=as.factor(sapply(strsplit(colnames(cutoff[,-c(1:6)]),'[0123456789]'),function(x) x[1])))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05)
# 176
cutoffc <- cutoff[p.value<0.05,]
sum((cutoffc$RT>=rt[21]&cutoffc$RT<=rt[30])|(cutoffc$RT>=rt[51]&cutoffc$RT<=rt[70]))
# 140
cutoffchange <- cutoff[cutoff$RT>=rt[21]&cutoff$RT<=rt[30]|cutoff$RT>=rt[51]&cutoff$RT<=rt[70],]
length(unique(paste(cutoffchange$mz,cutoffchange$rt)))
# 275
real <- read.csv('simcsv/control10.csv')
align <- enviGCMS::getalign(real$mz,cutoff$mz,real$rt,cutoff$RT)
length(unique(real$name[unique(align$xid)]))
# 99 
alignc <- enviGCMS::getalign(real$mz,cutoffc$mz,real$rt,cutoffc$RT)
length(unique(real$name[unique(alignc$xid)]))
# 29 detected changed compounds
alignchange <- enviGCMS::getalign(real$mz,cutoffchange$mz,real$rt,cutoffchange$RT)
length(unique(real$name[unique(alignchange$xid)]))
# 31 changed compounds
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 28 True positive 1 False positive 3 False negative
unique(real$name[unique(alignchange$xid)])



## -----------------------------------------------------------------------------------------------------
full <- read.csv('sim4mzmine.csv')
cutoff <- read.csv('simmzmine.csv')

fulldata <- full[,grepl('datafile(.*?)area',colnames(full))]
cutoffdata <- cutoff[,grepl('datafile(.*?)area',colnames(full))]
# impute NA just in case
fulldata[is.na(fulldata)] <- min(fulldata,na.rm = T)
cutoffdata[is.na(cutoffdata)] <- min(cutoffdata,na.rm = T)

library(genefilter)
rda <- rowttests(as.matrix(fulldata),fac=as.factor(sapply(strsplit(colnames(fulldata),'\\.|[0123456789]'),function(x) x[2])))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05,na.rm = T)
# 321
rt <- seq(10,590,length.out=100)

fullc <- full[p.value<0.05&(!is.na(p.value)),]
sum((fullc$rt*60>=rt[21]&fullc$rt*60<=rt[30])|(fullc$rt*60>=rt[51]&fullc$rt*60<=rt[70]))
# 302
fullchange <- full[full$rt*60>=rt[21]&full$rt*60<=rt[30]|full$rt*60>=rt[51]&full$rt*60<=rt[70],]
length(unique(paste(fullchange$mz,fullchange$rt)))
# 10422

# check compounds level
real <- read.csv('simcsv4/control10.csv')
align <- enviGCMS::getalign(real$mz,full$mz,real$rt,full$rt*60)
length(unique(real$name[unique(align$xid)]))
# 99
alignc <- enviGCMS::getalign(real$mz,fullc$mz,real$rt,fullc$rt*60)
length(unique(real$name[unique(alignc$xid)]))
# 25 detected changed compounds
alignchange <- enviGCMS::getalign(real$mz,fullchange$mz,real$rt,fullchange$rt*60)
length(unique(real$name[unique(alignchange$xid)]))
# 30 changed compounds
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 21 True positive 4 False positive 9 False negative

rda <- rowttests(as.matrix(cutoffdata),fac=as.factor(sapply(strsplit(colnames(cutoffdata),'\\.|[0123456789]'),function(x) x[2])))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05,na.rm = T)
# 53
rt <- seq(10,590,length.out=100)
cutoffc <- cutoff[p.value<0.05&(!is.na(p.value)),]
sum((cutoffc$rt*60>=rt[21]&cutoffc$rt*60<=rt[30])|(cutoffc$rt*60>=rt[51]&cutoffc$rt*60<=rt[70]))
# 51
cutoffchange <- cutoff[cutoff$rt*60>=rt[21]&cutoff$rt*60<=rt[30]|cutoff$rt*60>=rt[51]&cutoff$rt*60<=rt[70],]
length(unique(paste(cutoffchange$mz,fullchange$rt)))
# 10422

# check compounds level
real <- read.csv('simcsv4/control10.csv')
align <- enviGCMS::getalign(real$mz,cutoff$mz,real$rt,cutoff$rt*60)
length(unique(real$name[unique(align$xid)]))
# 99
alignc <- enviGCMS::getalign(real$mz,cutoffc$mz,real$rt,cutoffc$rt*60)
length(unique(real$name[unique(alignc$xid)]))
# 23 detected changed compounds
alignchange <- enviGCMS::getalign(real$mz,cutoffchange$mz,real$rt,cutoffchange$rt*60)
length(unique(real$name[unique(alignchange$xid)]))
# 30 changed compounds
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 22 True positive 1 False positive 8 False negative
