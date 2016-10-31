lv <- c(rep(1,4),rep(2,4),rep(3,4))
path <- "./data/acn/"
xset <- getopqedata(path)
raw <- svaplot(xset,lv)
praw <- svaplot(xset,lv,pqvalues = 'anova')
sv <- svaplot(xset,lv,pqvalues = 'sv')
psv <- svaplot(xset,lv,pqvalues = 'sv')

df <-  mutate(as.data.frame(sv),name = rownames(sv))
dreport <- annotateDiffreport(xset,metlin = T,polarity = "positive")

dreport$name <- paste0(round(dreport$mzmed,1),'/',round(dreport$rtmed))

a <- merge(df, dreport, by = "name",all.x = T)
ab <- a[(a$c18==4|a$c18==0)&(a$HLB==4|a$HLB==0)&(a$MM==4|a$MM==0),] 

c <- cbind(ab$rtmed,ab$mzmed)
plot(c)