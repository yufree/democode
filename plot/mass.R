source("http://bioconductor.org/biocLite.R")
biocLite("mzR")
library(mzR)
all <- openMSfile('./FULL200.CDF')
df <- header(all)
bb <- peaks(all)
aaaa <- sapply(bb,as.data.frame)
oddvals <- seq(1, ncol(aaaa), by=2) 
aaaaa <- unlist(aaaa[oddvals])
ccc <- unique(c(aaaaa))
ccc <- ccc[order(ccc)]
# bbb <- sapply(bb, "[",250:700)
# ddd <- unique(c(bbb))
# dddd <- ddd[ddd<700]
time <- df$retentionTime
df2 <- matrix(0, nrow = length(ccc), ncol = length(time))
rownames(df2) <- ccc
colnames(df2) <- time
rm(aaaa)
rm(aaaaa)
rm(oddvals)
rm(df)
rm(all)
gc()

for(i in 1:length(time)){
        temp <- bb[[i]]
        index <- which(ccc%in%temp[,1])
        df2[index,i] <- temp[,2]
}

ddd <- as.integer(ccc)
library(data.table)
dt = data.table(df2)
dt$fac <- ddd
df3 <- dt[,lapply(.SD, sum), by=ddd ]

df3 <- as.matrix(df3)
df7 <- df3[,2000:3000]
heatmap(df7)
library(rARPACK)
df4 <- svds(df3,2)
df5 <- df4$u %*% diag(df4$d) %*% t(df4$v)
rownames(df5) <- ddd
colnames(df5) <- time
df6 <- df5[,2000:3000]
heatmap(df6)
df8 <- as.data.frame(df5)
df9 <- as.data.frame(t(df8))
rownames(df8) <- ddd
colnames(df9) <- time
write.table(df3,'df3.txt')
