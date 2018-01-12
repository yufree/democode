# labeling report for control samples:
xcmsSet = xs3
sampleNames = sNctrl
intChoice = "intb"
groups = data.frame(xcmsSet@groups)
peakIntensities = getPeakIntensities(xcmsSet, 
                                     sampleNames, intChoice)
peakIntensities = peakIntensities[order(groups$mzmed), 
                                  ]
groups = groups[order(groups$mzmed), ]
groupRTs = groups$rtmed
groupMzs = groups$mzmed
groupIDs = as.numeric(rownames(groups))
nGroups = length(groupMzs)
classes = as.character(xcmsSet@phenoData[match(sampleNames, 
                                               rownames(xcmsSet@phenoData)), ])
numSamples = length(classes)
intensities1 = peakIntensities[, which(classes == 
                                               unlabeledSamples), drop = FALSE]
intensities2 = peakIntensities[, which(classes == 
                                               labeledSamples), drop = FALSE]
iMD = isotopeMassDiff
base = list()
labeled = list()
basePeak = list()
labeledPeak = list()
groupIndicesByRT = order(groupRTs)
orderedGroupRTs = groupRTs[groupIndicesByRT]
