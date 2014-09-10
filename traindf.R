# this function is used to divided a dataframe 'df' and response 'Y' into train 
# and test dataset for furthor modeling. N stand for the number of obs. in your 
# train data and if it is null, the default values is the 80 percent of the 
# whole obs.

traindf <- function (df,Y, N = 0)
{
        if(N==0){
                N <- round(nrow(df)*0.8)
        }
        train <- sample(nrow(df), N)
        trainX <- df[train, ]
        trainY <- Y[train]
        testX <- df[-train, ]
        testY <- Y[-train]
        return(list(trainX,trainY,testX,testY))
}