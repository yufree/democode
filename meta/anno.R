# devtools::install_github('yufree/MAIT')
# install.packages("WGCNA",repos="http://cran.r-project.org",dependencies=TRUE)
# install.packages("snow",repos="http://cran.r-project.org")
# install.packages("doSNOW",repos="http://cran.r-project.org")
# install.packages("parallel",repos="http://cran.r-project.org")
# install.packages("e1071",repos="http://cran.r-project.org")
# install.packages("XML",repos="http://cran.r-project.org")
# install.packages("R2HTML",repos="http://cran.r-project.org")
# install.packages("RCurl",repos="http://cran.r-project.org")
# source("http://bioconductor.org/biocLite.R")
# biocLite("Rdisop",suppressUpdates=TRUE)
# biocLite("SSOAP",suppressUpdates=TRUE)
# biocLite("KEGGREST",suppressUpdates=TRUE)
# install.packages("pcaMethods",repos="http://cran.r-project.org")
# install.packages("flashClust",repos="http://cran.r-project.org")
# install.packages("plyr",repos="http://cran.r-project.org")
# install.packages("png",repos="http://cran.r-project.org")
# install.packages("rjson",repos="http://cran.r-project.org")

source('https://raw.githubusercontent.com/yufree/democode/master/meta/getmetadata.R')
suppressPackageStartupMessages(library(MAIT))
suppressPackageStartupMessages(library(xMSannotator))
# use MAIT package
# Demo
# path <- "./data/"
# name <- "fishriver"
# anno(path,name)

anno <- function(path, name) {
        MAIT <- sampleProcessing(dataDir = path, project = name)
        MAIT <-
                peakAnnotation(
                        MAIT.object = MAIT,
                        corrWithSamp = 0.7,
                        corrBetSamp = 0.75,
                        perfwhm = 0.6
                )
        MAIT <-
                spectralSigFeatures(
                        MAIT.object = MAIT,
                        pvalue = 0.05,
                        p.adj = "BH",
                        scale = FALSE
                )
        signTable <-
                sigPeaksTable(MAIT.object = MAIT, printCSVfile = T)
        Biotransformations(MAIT.object = MAIT, peakPrecision = 0.005)
        MAIT <-
                identifyMetabolites(MAIT.object = MAIT, peakTolerance = 0.005)
        metTable <- metaboliteTable(MAIT)
        head(metTable)
        return(list(signTable,metTable))
}
# use xMSannotator package
# Demo
# path <- "./data/"
# xset <- getopqedata(path)
# result <- fanno(xset)
# other database options: KEGG, LipidMaps, T3DB
# queryadductlist=c("M+2H","M+H+NH4","M+ACN+2H","M+2ACN+2H","M+H","M+NH4","M+Na","M+ACN+H","M+ACN+Na","M+2ACN+H","2M+H","2M+Na","2M+ACN+H","M+2Na-H","M+H-H2O","M+H-2H2O") 
# other options: c("M-H","M-H2O-H","M+Na-2H","M+Cl","M+FA-H"); c("positive"); c("negative"); c("all");see data(adduct_table) for complete list
fanno <-
        function(xset,
                 outloc = "./result/",
                 mode = 'pos',
                 db_name = 'HMDB') {
                data(adduct_weights)
                data <- groupval(xset, 'medret', "into")
                mz <- groups(xset)[, 1]
                time <- groups(xset)[, 4]
                data <- as.data.frame(cbind(mz, time, data))
                data <- unique(data)
                annotres <-
                        multilevelannotation(
                                dataA = data,
                                max.mz.diff = 5,
                                max.rt.diff = 10,
                                cormethod = "pearson",
                                num_nodes = 12,
                                queryadductlist = c(
                                        "M+2H",
                                        "M+H+NH4",
                                        "M+ACN+2H",
                                        "M+2ACN+2H",
                                        "M+H",
                                        "M+NH4",
                                        "M+Na",
                                        "M+ACN+H",
                                        "M+ACN+Na",
                                        "M+2ACN+H",
                                        "2M+H",
                                        "2M+Na",
                                        "2M+ACN+H",
                                        "M+2Na-H",
                                        "M+H-H2O",
                                        "M+H-2H2O"
                                ),
                                mode = mode,
                                outloc = outloc,
                                db_name = db_name,
                                adduct_weights = adduct_weights,
                                num_sets = 1000,
                                allsteps = TRUE,
                                corthresh = 0.7,
                                NOPS_check = TRUE,
                                customIDs = NA,
                                missing.value = NA,
                                hclustmethod = "complete",
                                deepsplit = 2,
                                networktype = "unsigned",
                                minclustsize = 10,
                                module.merge.dissimilarity = 0.2,
                                filter.by = c("M+H"),
                                biofluid.location = NA,
                                origin = NA,
                                status = "Detected and Quantified",
                                boostIDs = NA,
                                max_isp = 5,
                                HMDBselect = "union",
                                mass_defect_window = 0.01,
                                pathwaycheckmode = "pm",
                                mass_defect_mode = mode
                        )
                return(annotres)
        }
mumdata <- function(xset, lv = NULL, name = 'test', method = "medret", intensity = 'inio'){
        data <- groupval(xset, method, intensity)
        if(intensity == "intb"){
                data[is.na(data)] = 0
        }
        if (is.null(lv)) {
                lv <- xset@phenoData[, 1]
        }
        mz <- xset@groups[, 1]
        rt <- xset@groups[, 4]
        mod <- model.matrix(~ lv)
        mod0 <- as.matrix(c(rep(1, ncol(data))))
        fstats <- sva::fstats(data,mod,mod0)
        pvalue <- sva::f.pvalue(data,mod,mod0)
        df <- cbind.data.frame(mz,rt,pvalue,fstats)
        filename <- paste0(name,'.txt')
        write.table(df, file = filename, sep = "\t",row.names = F)
        return(df)
}

# python2 main.py -c 0.05 -f test.txt -p 100 -m negative -o myoutput 
