## ----setup, include=FALSE-----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## -----------------------------------------------------------------------------------------------------
library(mzrtsim)
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
sn <- sample(c(100:10000),100)

# 30% changed
for(i in c(1:10)){
  simmzml(name=paste0('sim/case/case',i),db=monahrms1,pwidth = pw1,compound=compound,rtime = rt, sn=sn,unique = T)
}

for(i in c(1:10)){
  simmzml(name=paste0('sim/control/control',i),db=monahrms1,pwidth = pw2,compound=compound,rtime = rt, sn=sn,unique = T)
}


## -----------------------------------------------------------------------------------------------------
library(mzrtsim)
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
sn <- sample(c(100:10000),100)

# 30% changed
for(i in c(1:10)){
  simmzml(name=paste0('sim4/case/case',i),db=monahrms1,pwidth = pw1,compound=compound,rtime = rt, sn=sn,unique = T,inscutoff = 0)
}

for(i in c(1:10)){
  simmzml(name=paste0('sim4/control/control',i),db=monahrms1,pwidth = pw2,compound=compound,rtime = rt, sn=sn,unique = T,inscutoff = 0)
}

files <- list.files('sim4',full.names = T,recursive = T,pattern = '*.csv')
name <- basename(files)
file.rename(files,paste0('simcsv4/',name))


## -----------------------------------------------------------------------------------------------------
library(mzrtsim)
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
sn <- sample(c(100:10000),100)

for(i in c(1:10)){
  simmzml(name=paste0('sim5/tailing/tailing',i),db=monahrms1,pwidth = pw1,compound=compound,rtime = rt, sn=sn,unique = T,tailingfactor = 1.5,matrix = T)
}

for(i in c(1:10)){
  simmzml(name=paste0('sim5/normal/normal',i),db=monahrms1,pwidth = pw2,compound=compound,rtime = rt, sn=sn,unique = T,tailingfactor = 1,matrix = T)
}

for(i in c(1:10)){
  simmzml(name=paste0('sim5/leading/leading',i),db=monahrms1,pwidth = pw2,compound=compound,rtime = rt, sn=sn,unique = T,tailingfactor = 0.8,matrix = T)
}

files <- list.files('sim5',full.names = T,recursive = T,pattern = '*.csv')
name <- basename(files)
file.rename(files,paste0('simcsv5/',name))


## -----------------------------------------------------------------------------------------------------
library(xcms)
setMSnbaseFastLoad(opt = F)
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
raw_data <- readMSData(files = files,pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk")
cwp <- CentWaveParam(ppm=5,peakwidth = c(5,15))
xdata <- findChromPeaks(raw_data, param = cwp)
xdata <- groupChromPeaks(xdata, param = PeakDensityParam(pd$sample_group))

re <- enviGCMS::getmzrt(xdata,name='sim1')

compound <- read.csv('simcsv/case1.csv')
# 593 peaks
compoundsub <- compound[compound$mz>100&compound$mz<1000,]
# 533 peaks
align2 <- enviGCMS::getalign(compoundsub$mz,re$mz,compoundsub$rt,re$rt)
length(unique(compound$name[align2$xid]))
# 88
# 510 peaks
library(genefilter)
rda <- rowttests(as.matrix(re$data),fac=as.factor(re$group$sample_group))
p.value <- p.adjust(rda$p.value,'BH')

files <- list.files('sim',full.names = T,recursive = T,pattern = '*.csv')
name <- basename(files)
file.rename(files,paste0('simcsv/',name))


## -----------------------------------------------------------------------------------------------------
reticulate::use_python('/opt/homebrew/Caskroom/miniconda/base/bin/python')


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
data("monams1")
name <- sapply(monams1,function(x) x$name)
length(unique(name))
data("monahrms1")
name <- sapply(monahrms1,function(x) x$name)
length(unique(name))
data("hmdbcms")
name <- sapply(hmdbcms,function(x) x$name)
length(unique(name))
# peak number
pn <- sapply(monams1,function(x) x$np)
mean(as.numeric(pn))
median(as.numeric(pn))
pn <- sapply(hmdbcms,function(x) x$np)
mean(as.numeric(pn))
median(as.numeric(pn))


## -----------------------------------------------------------------------------------------------------
xcms <- enviGCMS::getmzrtcsv('sim5normalmzrt.csv')
openms <- read.csv('sim5openmsnormal.csv')
mzmine <- read.csv('sim5mzmingnormal.csv')

real <- read.csv('simcsv5/normal1.csv')

xcmsalign <- enviGCMS::getalign(real$mz,xcms$mz,real$rt,xcms$rt)
openmsalign <- enviGCMS::getalign(real$mz,openms$mz,real$rt,openms$RT)
mzminealign <- enviGCMS::getalign(real$mz,mzmine$mz,real$rt,mzmine$rt*60)

xcmsname <- paste(xcmsalign$mz2,xcmsalign$rt2)
openmsname <- paste(openmsalign$mz2,openmsalign$rt2)
mzminename <- paste(mzminealign$mz2,mzminealign$rt2)

length(unique(xcmsname))
length(unique(openmsname))
length(unique(mzminename))

library(ggvenn)
cvenn <- ggvenn(list(XCMS=real$name[unique(xcmsalign$xid)],OpenMS=real$name[unique(openmsalign$xid)],MZmine3=real$name[unique(mzminealign$xid)]))+ggtitle('B')
name <- paste(real$mz,real$rt)
pvenn <- ggvenn(list(XCMS=name[unique(xcmsalign$xid)],OpenMS=name[unique(openmsalign$xid)],MZmine3=name[unique(mzminealign$xid)]))+ggtitle('A')

library(patchwork)
p <- pvenn|cvenn
ggsave('figure1.png',p)

allfind <- c(real$name[unique(xcmsalign$xid)],real$name[unique(openmsalign$xid)],real$name[unique(mzminealign$xid)])
missing <- unique(real$name)[!unique(real$name)%in%allfind]

rm <- real[real$name%in%missing,]


## -----------------------------------------------------------------------------------------------------
xcms <- enviGCMS::getmzrtcsv('sim5leadingmzrt.csv')
openms <- read.csv('sim5openmsleading.csv')
mzmine <- read.csv('sim5mzmingleading.csv')

real <- read.csv('simcsv5/normal1.csv')

xcmsalign <- enviGCMS::getalign(real$mz,xcms$mz,real$rt,xcms$rt)
openmsalign <- enviGCMS::getalign(real$mz,openms$mz,real$rt,openms$RT)
mzminealign <- enviGCMS::getalign(real$mz,mzmine$mz,real$rt,mzmine$rt*60)

xcmsname <- paste(xcmsalign$mz2,xcmsalign$rt2)
openmsname <- paste(openmsalign$mz2,openmsalign$rt2)
mzminename <- paste(mzminealign$mz2,mzminealign$rt2)

length(unique(xcmsname))
length(unique(openmsname))
length(unique(mzminename))

ggvenn(list(XCMS=real$name[unique(xcmsalign$xid)],OpenMS=real$name[unique(openmsalign$xid)],MZmine3=real$name[unique(mzminealign$xid)]))+ggtitle('B')



## -----------------------------------------------------------------------------------------------------
xcms <- enviGCMS::getmzrtcsv('sim5tailingmzrt.csv')
openms <- read.csv('sim5openmstailing.csv')
mzmine <- read.csv('sim5mzminetailing.csv')

real <- read.csv('simcsv5/tailing1.csv')

xcmsalign <- enviGCMS::getalign(real$mz,xcms$mz,real$rt,xcms$rt)
openmsalign <- enviGCMS::getalign(real$mz,openms$mz,real$rt,openms$RT)
mzminealign <- enviGCMS::getalign(real$mz,mzmine$mz,real$rt,mzmine$rt*60)

xcmsname <- paste(xcmsalign$mz2,xcmsalign$rt2)
openmsname <- paste(openmsalign$mz2,openmsalign$rt2)
mzminename <- paste(mzminealign$mz2,mzminealign$rt2)

length(unique(xcmsname))
length(unique(openmsname))
length(unique(mzminename))

ggvenn(list(XCMS=real$name[unique(xcmsalign$xid)],OpenMS=real$name[unique(openmsalign$xid)],MZmine3=real$name[unique(mzminealign$xid)]))+ggtitle('B')



## -----------------------------------------------------------------------------------------------------
full <- enviGCMS::getmzrtcsv('sim4mzrt.csv')
cutoff <- enviGCMS::getmzrtcsv('sim1mzrt.csv')

library(genefilter)
rda <- rowttests(as.matrix(full$data),fac=as.factor(full$group$sample_group))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05)

rt <- seq(10,590,length.out=100)

fullc <- enviGCMS::getfilter(full,rowindex = p.value<0.05)
sum((fullc$rt>=rt[21]&fullc$rt<=rt[30])|(fullc$rt>=rt[51]&fullc$rt<=rt[70]))
# 817
fullchange <- enviGCMS::getfilter(full,rowindex = full$rt>=rt[21]&full$rt<=rt[30]|full$rt>=rt[51]&full$rt<=rt[70])

real <- read.csv('simcsv4/case1.csv')
align <- enviGCMS::getalign(real$mz,full$mz,real$rt,full$rt)
length(unique(real$name[unique(align$xid)]))
# 94
alignc <- enviGCMS::getalign(real$mz,fullc$mz,real$rt,fullc$rt)
length(unique(real$name[unique(alignc$xid)]))
# 41
alignchange <- enviGCMS::getalign(real$mz,fullchange$mz,real$rt,fullchange$rt)
length(unique(paste(fullchange$mz,fullchange$rt)))
# 914
length(unique(real$name[unique(alignchange$xid)]))
# 26
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 41-26 = 15

rda <- rowttests(as.matrix(cutoff$data),fac=as.factor(cutoff$group$sample_group))
p.value <- p.adjust(rda$p.value,'BH')
sum(p.value<0.05)
# 198
cutoffc <- enviGCMS::getfilter(cutoff,rowindex = p.value<0.05)
sum((cutoffc$rt>=rt[21]&cutoffc$rt<=rt[30])|(cutoffc$rt>=rt[51]&cutoffc$rt<=rt[70]))
# 154
cutoffchange <- enviGCMS::getfilter(cutoff,rowindex = (cutoff$rt>=rt[21]&cutoff$rt<=rt[30])|(cutoff$rt>=rt[51]&cutoff$rt<=rt[70]))

real <- read.csv('simcsv/case1.csv')
align <- enviGCMS::getalign(real$mz,cutoff$mz,real$rt,cutoff$rt)
length(unique(real$name[unique(align$xid)]))
# 98
alignc <- enviGCMS::getalign(real$mz,cutoffc$mz,real$rt,cutoffc$rt)
length(unique(real$name[unique(alignc$xid)]))
# 35
alignchange <- enviGCMS::getalign(real$mz,cutoffchange$mz,real$rt,cutoffchange$rt)
length(unique(paste(cutoffchange$mz,cutoffchange$rt)))
# 154
length(unique(real$name[unique(alignchange$xid)]))
# 27
sum(unique(real$name[unique(alignc$xid)]) %in% unique(real$name[unique(alignchange$xid)]))
# 35-27 = 8

