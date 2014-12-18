#fit regression and plot to show that it works
load('data.RData')
reg <- lm(Y~X,data=data)
plot(reg$fitted,reg$resid,pch=16,main="Residual plot from data")
