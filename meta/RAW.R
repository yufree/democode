getRAW <- function(xset,lv){
        lsa <- list()
        index <- 0
        m <- length(table(lv))
        for(i in 1:m){
                if(i == 1){
                        index <- index
                }else{
                        index <- index+table(lv)[i-1]
                }
                ls <- list()
                for(j in 1:table(lv)[i]){
                        RAW <- getXcmsRaw(xset,sampleidx=j+index)
                        RAWdata <- RAW@env$profile
                        colnames(RAWdata) <- RAW@scantime
                        rownames(RAWdata) <- seq(RAW@mzrange[1],RAW@mzrange[2],0.1)
                        ls[[j]] <- RAWdata
                }
                lsa[[i]] <- ls
        }
        return(lsa)
}

getalign <- function(list){
        a <- lapply(list,rownames)
        c <- Reduce(intersect,a)
        e <- lapply(list,function(x){x[c,]})
        return(e)
}
