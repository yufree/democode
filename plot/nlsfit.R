library(minpack.lm)
x <- c(0,2,4,6,8,10,12)
y1 <- c(0,304,2101,2003,1622,2056,1584)
y2 <- c(0,351,187,2674,4236,4673,3190,2796)
p1 <- 1500
p2 <- 0.5
fit = nls(y1 ~ p1*(1-exp(-p2*x)), start=list(p1=p1,p2=p2))
fit = nlsLM(y1 ~ p1*(1-exp(-p2*x)), start=list(p1=p1,p2=p2))
x0 <- 1:100
y <- 1941*(1-exp(-0.36*x))
plot(y1~x)
lines(y~x)
