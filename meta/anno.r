suppressWarnings(suppressPackageStartupMessages(library(xcms)))
library(RColorBrewer)
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(limma))
library(CAMERA)
library(qvalue)
library(MAIT)

a <- MAITbuilder(
        data = df$dataCorrected,
        spectraID = NULL,
        masses = df$mz,
        rt = df$rt,
        spectraEstimation = TRUE,
        classes = lv,
        rtRange = 0.2,
        corThresh = 0.7
)
importMAIT <- Biotransformations(
        MAIT.object = a,
        adductAnnotation = TRUE,
        peakPrecision = 0.005,
        adductTable = NULL
)
importMAIT <- identifyMetabolites(MAIT.object = importMAIT,
                                  peakTolerance = 0.005,
                                  polarity = "positive")
