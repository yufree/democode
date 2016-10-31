data1 <- read.csv('data/cn.csv')
data2 <- read.csv('data/trophic.csv')
png('data1.png',width = 8, height = 8, units = 'in', res = 300)
par(mar = c(5,5,4,2))
plot(data1$N~data1$C,
     xlim = c(-40,-10),ylim = c(-2,16),
     pch = as.numeric(data1$Cat),
     col = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8),
     xlab = expression(paste(delta^{13}, "C (\u2030)")),ylab = expression(paste(delta^{15}, "N (\u2030)")))

legend('topleft',legend = c(expression(italic('O. latipes')),expression(italic('C. japonica')),expression(italic('P. acuta')),expression(italic('C. fluminea')),expression(italic('D. magna')),expression(italic('P. stratiotes')),expression(italic('S. borealis')),"Sediment"),pch = unique(as.numeric(data1$Cat)),col = c(1:8))
dev.off()


# Make up some data
new <- data.frame(tro = seq(0, 5, length.out = 200))

# fit a quadratic
model <- lm(con~tro,data = data2)
fitted <- predict(model,newdata = new, interval = "prediction")
png('data2.png',width = 8, height = 8, units = 'in', res = 300)
# plot the data and the fitted line
par(mar = c(5,5,4,2))
plot(data2$con~data2$tro,
     pch = as.numeric(data2$specious),
     col = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7),
     xlim=c(0,5),ylim = c(0,10),
     xlab = 'Trophic level',
     ylab = expression(paste("Log"[Ceerium],'(',mu,'g/kg), lipid')))
legend('topright',legend = c(expression(italic('O. latipes')),expression(italic('C. japonica')),expression(italic('P. acuta')),expression(italic('C. fluminea')),expression(italic('D. magna')),expression(italic('P. stratiotes')),expression(italic('S. borealis'))),pch = unique(as.numeric(data1$Cat)),col = c(1:8))


abline(model,lwd=2)

matlines(new$tro,fitted[,c("lwr","upr")],col=2,lty=2,type="l",pch=19,lwd = 2)

text(1,1,'R = -0.5940, TMF = 0.6527',cex=1.2)
dev.off()
