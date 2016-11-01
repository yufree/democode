con <- c(1.202,1.221,1.241,1.26,1.279,1.297,1.316)

confix <- c(1.301,1.38,1.398,1.398,rep(1.447,5),1.519,1.531)
tl <- c(2.03,1.48,1.15,1.25,1.37,1.83,1.52,2.46,1.63,3.2,2.93,3.01,3.16,3.1,2.19,1.91,2.99,2.71)

mean <- mean(con)
sd <- sd(con)

int <- NULL
slp <- NULL

for(i in 1:10000){
        scon <- rnorm(7,mean =mean,sd =sd)
        con0 <- c(scon,confix)
        fit <- summary(lm(tl~con0))
        lmint <- fit$coefficients[1,1]
        lmslp <- fit$coefficients[2,1]
        int <- c(int,lmint)
        slp <- c(slp,lmslp)
}

mean(int)
mean(slp)
sd(int)
sd(slp)

t.test(int)
t.test(slp)
