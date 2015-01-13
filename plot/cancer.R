cancer <- read.csv('data/cancer.csv')

fit0 <- lm(log(cancer$Lifetime.cancer.incidence)~log(cancer$Cumulative.number.of.divisions.of.all.stem.cells.per.lifetime..lscd.))

fit1 <- lm(cancer$Lifetime.cancer.incidence~cancer$Cumulative.number.of.divisions.of.all.stem.cells.per.lifetime..lscd.)

fit2 <- lm(cancer$Lifetime.cancer.incidence~log(cancer$Cumulative.number.of.divisions.of.all.stem.cells.per.lifetime..lscd.))

fit3 <- lm(cancer$Lifetime.cancer.incidence~sqrt(cancer$Cumulative.number.of.divisions.of.all.stem.cells.per.lifetime..lscd.))   

summary(fit0)$r.squared
summary(fit1)$r.squared
summary(fit2)$r.squared
summary(fit3)$r.squared

png('cancer.png')
par(mfrow=c(2,2))

plot(log(cancer$Lifetime.cancer.incidence)~log(cancer$Cumulative.number.of.divisions.of.all.stem.cells.per.lifetime..lscd.),pch=19,main=expression(paste('log x and log y ',r^2,' 0.646')),xlab='Total stem cell division',ylab='lifetime risk')

plot(cancer$Lifetime.cancer.incidence~cancer$Cumulative.number.of.divisions.of.all.stem.cells.per.lifetime..lscd.,pch=19,main=expression(paste('no transformation ',r^2,' 0.284')),xlab='Total stem cell division',ylab='lifetime risk')

plot(cancer$Lifetime.cancer.incidence~log(cancer$Cumulative.number.of.divisions.of.all.stem.cells.per.lifetime..lscd.),pch=19,main=expression(paste('log x ',r^2,' 0.183')),xlab='Total stem cell division',ylab='lifetime risk')

plot(cancer$Lifetime.cancer.incidence~sqrt(cancer$Cumulative.number.of.divisions.of.all.stem.cells.per.lifetime..lscd.),pch=19,main=expression(paste('squared x ',r^2,' 0.356')),xlab='Total stem cell division',ylab='lifetime risk')
dev.off()

