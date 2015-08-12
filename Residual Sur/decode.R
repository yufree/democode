# this codes partly came from Prof. John Staudenmayer's code which is post on Prof. Leonard A. Stefanski's website.
# fit regression and plot to show that it works
load('data.RData')
reg <- lm(X1~ .,data=data)
plot(reg$fitted,reg$resid,pch=16,main="Residual plot from data")
