set.seed(42)
signal1 <- rep(75,5)
batch1 <- c(rep(30,3),rep(-20,2))
error1 <- rnorm(5,sd = 10)
real1 <- signal1 + batch1 + error1
cor1 <- signal1 + error1

signal2 <- rep(150,5)
batch2 <- c(rep(-20,3),rep(50,2))
error2 <- rnorm(5,sd = 10)
real2 <- signal2 + batch2 + error2
cor2 <- signal2 + error2

signal <- c(signal1,signal2)

batch <- c(batch1,batch2)

error <- c(error1,error2)

data <- c(real1,real2)

cdata <- c(cor1,cor2)

lv <- factor(c(rep('1',5),rep('2',5)))

mod <- stats::model.matrix(~lv)
fit <- limma::lmFit(t(data),mod)
signal0 <- fit$coef[, 1:nlevels(lv)] %*% t(mod[,1:nlevels(lv)])
signal0 <- as.numeric(t(signal0))
error0 <- data - signal0
data0 <- cbind(data,signal0,error0)
colnames(data0) <- c('raw data','signal','error')

data <- cbind(data,signal,batch,error,cdata)
colnames(data) <- c('raw data','signal','batch','error','corrected data')

pdf('scheme2.pdf',width = 9,height = 7)
# par(mar = c(3,0,3,0),mfrow = c(2,1))
par(mar = c(3,0,3,0))
# barplot(data0,beside = T,col = c(rep(rgb(0,0,0,0.5),5),rep(rgb(0,0,1,0.5),5)),yaxt = 'n',main = 'without correction')
# legend('topright', legend = c('Group 1', 'Group 2'), pch = c(15,15),bty = 'n', col = c(rgb(0,0,0,0.5),rgb(0,0,1,0.5)))

barplot(data,beside = T,col = c(rep(rgb(0,0,0,0.5),5),rep(rgb(0,0,1,0.5),5)),yaxt = 'n')
legend('topright', legend = c('Group 1', 'Group 2'), pch = c(15,15),bty = 'n', col = c(rgb(0,0,0,0.5),rgb(0,0,1,0.5)))
dev.off()
