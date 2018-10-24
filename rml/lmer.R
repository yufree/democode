# build the model
anovamixed <- function(x){
    ls1 <- lme4::lmer(x~time+(time|subject))
    ls2 <- lme4::lmer(x~1+(time|subject))
    # ls1 <- lme4::lmer(x~time+(time|subject))
    # ls2 <- lme4::lmer(x~1+(time|subject))
    ls <- summary(ls1)
    t <- anova(ls1,ls2)
    return(c(t$`Pr(>Chisq)`[2],ls$coefficients[2,1]))
}
a <- apply(dataf,1,anovamixed)
# FDR control
padjust <- p.adjust(a[1,],method = 'BH')

sum(padjust<0.05&a[2,]>0)
sum(padjust<0.05&a[2,]<0)
# get the peaks with changes
datap <- dataf[padjust<0.05&a[2,]>0,]
datan <- dataf[padjust<0.05&a[2,]<0,]
rownames(datan) <- rownames(dataf)[padjust<0.05&a[2,]<0]
# plot the compounds
i = 0
apply(datan,1,visindividual)
i = 0
apply(datan,1,visindividual2)

i = 0
apply(datan,1,visindividualall)

# for multiple group nlme model and table output

result <- ndf %>%
        as.tibble() %>%
        nest(-group) %>%
        mutate(model = map(data,~ lme(values ~ factor1+factor2+factor3-1,random=~1|sample,data = .)),tidy = map(model,tidy,effects = 'fixed')) %>%
        unnest(tidy)

result <- ndf %>%
        as.tibble() %>%
        nest(-group) %>%
        mutate(model = map(data,~ lme(values ~ factor1+factor2+factor3-1,random=~1|sample,data = .)))

names(result$model) <- result$group

table <- result$model %>%
        huxreg(tidy_args = list(effects = 'fixed'))

quick_docx(table,file = 'table.docx')