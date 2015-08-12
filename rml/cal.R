data <- read.csv('data/model.csv',header = F)
datatest <- read.csv('data/modeltest.csv',header = F)
N = nrow(data)
k = 10
folds = sample(1:k, N, replace = TRUE)

cv.errors = NA
listmodel <- NA

for (j in 1:k) {
        best.fit = lm(data[folds != j,4]~data[folds != j,2])
        pred = predict(best.fit,data.frame(data[folds == j,2]))
        listmodel[j] = best.fit
        cv.errors[j] = mean((data[folds == j,4]-pred)^2)
}

minX <- which.min(cv.errors)
md1 <- listmodel[minX]
min(cv.errors)
# test
predtest <- datatest[,2]*md1[[1]][2]+md1[[1]][1]
cv.errorstest <- mean((datatest[,4]-predtest)^2,na.rm = T)
cv.errorstest
# compare
fit <- lm(data[,4]~data[,2])
pred0 = datatest[,2]*coefficients(fit)[2]+coefficients(fit)[1]
error0 = mean((data[,4]-data[,2]*coefficients(fit)[2]-coefficients(fit)[1])^2)
error0
error = mean((datatest[,4]-pred0)^2,na.rm = T)
error

cv.errors = NA
listmodel <- NA

for (j in 1:k) {
        best.fit = lm(data[folds != j,4]~data[folds != j,3])
        pred = predict(best.fit,data.frame(data[folds == j,3]))
        listmodel[j] = best.fit
        cv.errors[j] = mean((data[folds == j,4]-pred)^2)
}

minX <- which.min(cv.errors)
md2 <- listmodel[minX]
min(cv.errors)
# test
predtest <- datatest[,3]*md2[[1]][2]+md2[[1]][1]
cv.errorstest <- mean((datatest[,4]-predtest)^2,na.rm = T)
cv.errorstest
# compare
fit <- lm(data[,4]~data[,3])
pred0 = datatest[,3]*coefficients(fit)[2]+coefficients(fit)[1]
error0 = mean((data[,4]-data[,3]*coefficients(fit)[2]-coefficients(fit)[1])^2)
error0
error = mean((datatest[,4]-pred0)^2,na.rm = T)
error

cv.errors = NA
listmodel <- NA

for (j in 1:k) {
        best.fit = lm(data[folds != j,4]~data[folds != j,5])
        pred = predict(best.fit,data.frame(data[folds == j,5]))
        listmodel[j] = best.fit
        cv.errors[j] = mean((data[folds == j,4]-pred)^2)
}

minX <- which.min(cv.errors)
md3 <- listmodel[minX]
min(cv.errors)
# test
predtest <- datatest[,5]*md3[[1]][2]+md3[[1]][1]
cv.errorstest <- mean((datatest[,4]-predtest)^2,na.rm = T)
cv.errorstest
# compare
fit <- lm(data[,4]~data[,5])
pred0 = datatest[,5]*coefficients(fit)[2]+coefficients(fit)[1]
error0 = mean((data[,4]-data[,5]*coefficients(fit)[2]-coefficients(fit)[1])^2)
error0
error = mean((datatest[,4]-pred0)^2,na.rm = T)
error

predtest0 <- (datatest[,2]*md1[[1]][2]+md1[[1]][1]+datatest[,3]*md2[[1]][2]+md2[[1]][1]+datatest[,5]*md3[[1]][2]+md3[[1]][1])/3
errormix <- mean((datatest[,4]-predtest0)^2,na.rm = T)
errormix

predgq <- datatest[,2]*1.004139409+0.2060091161
predgq2 <- data[,2]*1.30929533+0.08901246794
mean((datatest[,4]-predgq)^2,na.rm = T)

plot(datatest[,1],datatest[,4],type='l')
points(datatest[,1],predtest0,col='red',type='l')
t.test(datatest[1:3300,4],predtest0[1:3300],paired = T)
#points(datatest[,1],datatest[,5]*md3[[1]][2]+md3[[1]][1],col='blue',type='l')
points(datatest[,1],predgq,col='green',type='l')
points(datatest[,1],datatest[,2]*md1[[1]][2]+md1[[1]][1],col='blue',type='l')
points(datatest[,1],datatest[,3]*md2[[1]][2]+md2[[1]][1],col='yellow',type='l')
points(data[,1],data[,5]*md3[[1]][2]+md3[[1]][1],col='red',type='l')
points(data[,1],predgq2,col='purple',type='l')
points(datatest[,1],pred0,col='red',type='l')
