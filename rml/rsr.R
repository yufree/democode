df0 <- read.csv('./2008_BOLD_Survey_of_Puget_Sound_Levels_for_PCB_Congeners_Raw_Data.csv')
df <- df0[,8:173]
rownames(df) <- df0$Station.ID
dfm <- data.matrix(df)
dfm <- round(dfm)
idx <- sample(1:length(dfmlog),3000)
dfmna[idx] <- NA
nam <- is.na(dfmna)
nam <- 1-as.numeric(nam)
nam <- matrix(nam,nrow = 75)
library(R.matlab)
writeMat('dfmo.mat', df = dfm, Y = dfmna, R = nam)

hist(dfmlog <- log(dfm))
image(dfmlog, main = "Raw data")
set.seed(42)
idx <- sample(1:length(dfmlog),3000)
dfmlogna <- dfmlog
dfmlogna[idx] <- NA
image(dfmlog)
image(dfmlogna)
dfmna <- dfm 
dfmna[idx] <- NA
nam <- is.na(dfmna)
nam <- 1-as.numeric(nam)
l  <- cut(dfmlog, c(-Inf, -5, -3, 0, 5, 10), include.lowest=TRUE, labels = FALSE)
dft <- matrix(l,nrow = 75)
rownames(dft) <- df0$Station.ID
colnames(dft) <- colnames(df)
idx <- sample(1:length(dft),3000)
dftna <- dft 
dftna[idx] <- 0
image(dft)
image(dftna)

library(recommenderlab)
dfr <- as(dfmna, "realRatingMatrix")
# dfr <- as(t(dfmlogna), "realRatingMatrix")
scheme <- evaluationScheme(dfr,method = "split", train = .9,given = 10, goodRating = 1)
algorithms <- list(
        "random items" = list(name="RANDOM", param=list(normalize = "Z-score")),
        "popular items" = list(name="POPULAR", param=list(normalize = "Z-score")),
        "user-based CF" = list(name="UBCF", param=list(normalize = "Z-score",
                                                       method="Cosine",
                                                       nn=50, minRating=3)),
        "item-based CF" = list(name="IBCF2", param=list(normalize = "Z-score"
        ))
        
)
results <- evaluate(scheme, algorithms, n=c(1, 3, 5, 10, 15, 20))

# Draw ROC curve
plot(results, annotate = 1:4, legend="topleft")

# See precision / recall
plot(results, "prec/rec", annotate=3)

r <- Recommender(dfr, method="UBCF", param=list(normalize = "Z-score",
                                                method="Cosine",
                                                nn=50, minRating=3))

recom <- predict(r, dfr, type = 'ratings')

as(recom, "list")

dfmcv <- as(recom, "matrix")
dftt <- dfm
dftt[-idx] <- NA
dfmerror <- sum(abs(dfmcv - dftt),na.rm = T)
image(dfmcv)
image(dftt)
recom3 <- bestN(recom,n=3)

as(recom3, "list")
library(R.matlab)

writeMat('dfmo.mat', dfm = dfm)
writeMat('dfm.mat', dfm = dft, Y = dftna)

