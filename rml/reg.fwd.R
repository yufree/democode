library(leaps)

# N stand for the maximum size of subsets to examine
# this is a 5-fold cross-validation

source('traindf.R')
source('predict.regsubsets.R')

data <- traindf(df,Y)

trainX <- data[[1]]
trainY <- data[[2]]
testX <- data[[3]]
testY <- data[[4]]

N = 50

k = 5

folds = sample(1:k, nrow(traindataX), replace = TRUE)

cv.errors = matrix(NA, k, N, dimnames = list(NULL, paste(1:N)))

for (j in 1:k) {
        best.fit = regsubsets(y = trainY[folds != j],
                              x = trainX[folds != j, ],
                              nvmax = N, 
                              really.big = T,
                              method = "forward"
                              )
        for (i in 1:50) {
                pred = predict.regsubsets(best.fit,
                                          traindX[folds == j, ],
                                          id = i)
                cv.errors[j, i] = mean((trainY[folds == j] - pred)^2)
        }
}

mean.cv.errors = apply(cv.errors, 2, mean)

minX <- which.min(mean.cv.errors)

reg.fwd = regsubsets(x = trainX, 
                     y = trainY, 
                     nvmax = minX, 
                     really.big = T, 
                     method = "forward")

val.errors <- rep(NA, minX)

for (i in 1:minX) {
        pred = predict.regsubsets(reg.fwd, 
                                  testX, 
                                  id = i)
        val.errors[i] = mean((testY - pred)^2)}

minX <- which.min(val.errors)

reg.fwd.final = regsubsets(x = df, 
                     y = Y, 
                     nvmax = minX, 
                     really.big = T, 
                     method = "forward")

coef(reg.fwd.final, minX)
