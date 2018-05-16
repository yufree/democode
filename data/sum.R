t <- c(10, 30, 50, 70, 100)
r <- c(1, 5, 10, 25, 30, 50)
m <- c(1:10)
re2 <- NULL
for (i in c(1:length(m))) {
        re <- val <- NULL
        for (k in c(1:length(r))) {
                for (j in c(1:length(t))) {
                        result <- NULL
                        for (l in c(1:m[i])) {
                                a <- 1 / (l ^ 2) * exp((-1020 * 9 * pi ^ 2 * l^2 * t[j]) / (16 * r[k] ^ 2))
                                result <- c(result, a)
                        }
                        print(sum(result))
                        val <- 1 - 64 / (9 * pi ^ 2) * sum(result)
                        re <- rbind(re, cbind(val,t[j],r[k]))
                }
        }
        re2 <- rbind(re2,cbind(re,m[i]))
}
colnames(re2) <- c('D','t','r','n')
write.csv(re2,'resulteq.csv')
