# devtools::install_github('yufree/sva-devel')
suppressWarnings(suppressPackageStartupMessages(library(xcms)))
library(RColorBrewer)
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(MAIT))
library(CAMERA)
library(qvalue)

svacor <-
        function(xset,
                 lv = NULL,
                 annotation = F,
                 polarity = "positive",
                 nSlaves = 12) {
                data <- groupval(xset, "maxint", value = 'into')
                if (is.null(lv)) {
                        lv <- xset@phenoData[, 1]
                }
                mz <- xset@groups[, 1]
                rt <- xset@groups[, 4]
                mod <- model.matrix( ~ lv)
                mod0 <- as.matrix(c(rep(1, ncol(data))))
                svafit <- sva(data, mod)
                if (svafit$n.sv == 0) {
                        svaX <- model.matrix( ~ lv)
                        lmfit <- lmFit(data, svaX)
                        signal <-
                                lmfit$coef[, 1:nlevels(lv)] %*% t(svaX[, 1:nlevels(lv)])
                        error <- data - signal
                        rownames(signal) <-
                                rownames(error) <- rownames(data)
                        colnames(signal) <-
                                colnames(error) <- colnames(data)
                        pValues = f.pvalue(data, mod, mod0)
                        qValues = qvalue(pValues)
                        qValues = qValues$qvalues
                        if (annotation) {
                                dreport <-
                                        annotateDiffreport(
                                                xset,
                                                metlin = T,
                                                polarity = polarity,
                                                nSlaves = nSlaves
                                        )
                                dreport <-
                                        dreport[order(as.numeric(rownames(dreport))), ]
                                li <-
                                        list(data,
                                             signal,
                                             error,
                                             pValues,
                                             qValues,
                                             dreport,
                                             mz,
                                             rt)
                                names(li) <-
                                        c(
                                                'data',
                                                'signal',
                                                'error',
                                                'p-values',
                                                'q-values',
                                                'diffreport',
                                                'mz',
                                                'rt'
                                        )
                        } else{
                                li <- list(data,
                                           signal,
                                           error,
                                           pValues,
                                           qValues,
                                           mz,
                                           rt)
                                names(li) <-
                                        c(
                                                'data',
                                                'signal',
                                                'error',
                                                'p-values',
                                                'q-values',
                                                'mz',
                                                'rt'
                                        )
                        }
                }
                else{
                        message('Data is correcting ...')
                        svaX <- model.matrix( ~ lv + svafit$sv)
                        lmfit <- lmFit(data, svaX)
                        batch <-
                                lmfit$coef[, (nlevels(lv) + 1):(nlevels(lv) + svafit$n.sv)] %*% t(svaX[, (nlevels(lv) +
                                                                                                                  1):(nlevels(lv) + svafit$n.sv)])
                        signal <-
                                lmfit$coef[, 1:nlevels(lv)] %*% t(svaX[, 1:nlevels(lv)])
                        error <- data - signal - batch
                        datacor <- signal + error
                        svaX2 <- model.matrix( ~ lv)
                        lmfit2 <- lmFit(data, svaX2)
                        signal2 <-
                                lmfit2$coef[, 1:nlevels(lv)] %*% t(svaX2[, 1:nlevels(lv)])
                        error2 <- data - signal2
                        rownames(signal2) <-
                                rownames(error2) <-
                                rownames(datacor) <-
                                rownames(signal) <-
                                rownames(batch) <-
                                rownames(error) <- rownames(data)
                        colnames(signal2) <-
                                colnames(error2) <-
                                colnames(datacor) <-
                                colnames(signal) <-
                                colnames(batch) <-
                                colnames(error) <- colnames(data)
                        
                        modSv = cbind(mod, svafit$sv)
                        mod0Sv = cbind(mod0, svafit$sv)
                        pValuesSv = f.pvalue(data, modSv, mod0Sv)
                        qValuesSv = qvalue(pValuesSv)
                        qValuesSv = qValuesSv$qvalues
                        
                        pValues = f.pvalue(data, mod, mod0)
                        qValues = qvalue(pValues)
                        qValues = qValues$qvalues
                        if (annotation) {
                                dreport <-
                                        annotateDiffreport(
                                                xset,
                                                metlin = T,
                                                polarity = polarity,
                                                nSlaves = nSlaves
                                        )
                                dreport <-
                                        dreport[order(as.numeric(rownames(dreport))), ]
                                li <-
                                        list(
                                                data,
                                                datacor,
                                                signal,
                                                batch,
                                                error,
                                                signal2,
                                                error2,
                                                pValues,
                                                qValues,
                                                pValuesSv,
                                                qValuesSv,
                                                dreport,
                                                svafit$pprob.gam,
                                                svafit$pprob.b,
                                                mz,
                                                rt
                                        )
                                names(li) <-
                                        c(
                                                'data',
                                                'dataCorrected',
                                                'signal',
                                                'batch',
                                                'error',
                                                'signal2',
                                                'error2',
                                                'p-values',
                                                'q-values',
                                                'p-valuesCorrected',
                                                'q-valuesCorrected',
                                                'diffreport',
                                                'PosteriorProbabilitiesSurrogate',
                                                'PosteriorProbabilitiesMod',
                                                'mz',
                                                'rt'
                                        )
                        }
                        else{
                                li <-
                                        list(
                                                data,
                                                datacor,
                                                signal,
                                                batch,
                                                error,
                                                signal2,
                                                error2,
                                                pValues,
                                                qValues,
                                                pValuesSv,
                                                qValuesSv,
                                                svafit$pprob.gam,
                                                svafit$pprob.b,
                                                mz,
                                                rt
                                        )
                                names(li) <-
                                        c(
                                                'data',
                                                'dataCorrected',
                                                'signal',
                                                'batch',
                                                'error',
                                                'signal2',
                                                'error2',
                                                'p-values',
                                                'q-values',
                                                'p-valuesCorrected',
                                                'q-valuesCorrected',
                                                'PosteriorProbabilitiesSurrogate',
                                                'PosteriorProbabilitiesMod',
                                                'mz',
                                                'rt'
                                        )
                        }
                        message('Done!')
                }
                return(li)
        }

svapca <- function(list,
                   center = T,
                   scale = T) {
        data <- list$data
        Signal <- list$signal
        Batch <- list$batch
        error <- list$error
        datacor <- list$dataCorrected
        
        par(mfrow = c(2, 5), mar = c(4, 4, 2.6, 1))
        
        pcao <- prcomp(t(data), center = center, scale = scale)
        pcaoVars = signif(((pcao$sdev) ^ 2) / (sum((pcao$sdev) ^ 2)), 3) *
                100
        plot(pcao, type = "l", main = "PCA")
        
        pca <- prcomp(t(Signal), center = TRUE, scale = TRUE)
        pcaVars = signif(((pca$sdev) ^ 2) / (sum((pca$sdev) ^ 2)), 3) *
                100
        plot(pca, type = "l", main = "PCA-signal")
        
        pcab <- prcomp(t(Batch), center = center, scale = scale)
        pcabVars = signif(((pcab$sdev) ^ 2) / (sum((pcab$sdev) ^ 2)), 3) *
                100
        plot(pcab, type = "l", main = "PCA-batch")
        
        pcae <- prcomp(t(datacor), center = center, scale = scale)
        pcaeVars = signif(((pcae$sdev) ^ 2) / (sum((pcae$sdev) ^ 2)), 3) *
                100
        plot(pcae, type = "l", main = "PCA-error")
        
        pcac <- prcomp(t(datacor), center = center, scale = scale)
        pcacVars = signif(((pcac$sdev) ^ 2) / (sum((pcac$sdev) ^ 2)), 3) *
                100
        plot(pcac, type = "l", main = "PCA-corrected")
        
        plot(
                pcao$x[, 1],
                pcao$x[, 2],
                xlab = paste("PC1:", pcaoVars[1], "% of Variance Explained"),
                ylab = paste("PC2:", pcaoVars[2], "% of Variance Explained"),
                pch = colnames(data),
                cex = 2,
                main = "PCA"
        )
        
        plot(
                pca$x[, 1],
                pca$x[, 2],
                xlab = paste("PC1:", pcaVars[1], "% of Variance Explained"),
                ylab = paste("PC2:", pcaVars[2], "% of Variance Explained"),
                pch = colnames(Signal),
                cex = 2,
                main = "PCA-signal"
        )
        
        plot(
                pcab$x[, 1],
                pcab$x[, 2],
                xlab = paste("PC1:", pcabVars[1], "% of Variance Explained"),
                ylab = paste("PC2:", pcabVars[2], "% of Variance Explained"),
                pch = colnames(Batch),
                cex = 2,
                main = "PCA-batch"
        )
        
        plot(
                pcae$x[, 1],
                pcae$x[, 2],
                xlab = paste("PC1:", pcaeVars[1], "% of Variance Explained"),
                ylab = paste("PC2:", pcaeVars[2], "% of Variance Explained"),
                pch = colnames(error),
                cex = 2,
                main = "PCA-error"
        )
        
        plot(
                pcac$x[, 1],
                pcac$x[, 2],
                xlab = paste("PC1:", pcacVars[1], "% of Variance Explained"),
                ylab = paste("PC2:", pcacVars[2], "% of Variance Explained"),
                pch = colnames(datacor),
                cex = 2,
                main = "PCA-corrected"
        )
}

svaplot <- function(list,
                    pqvalues = "sv",
                    pt = 0.05,
                    qt = 0.05) {
        data <- list$data
        signal <- list$signal
        signal2 <- list$signal2
        batch <- list$batch
        error <- list$error
        error2 <- list$error2
        datacor <- list$dataCorrected
        pValues <- list$'p-values'
        qValues <- list$'q-values'
        pValuesSv <- list$'p-valuesCorrected'
        qValuesSv <- list$'q-valuesCorrected'
        
        icolors <-
                colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)
        
        if (is.null(signal2)) {
                if (pqvalues == "anova" & sum(pValues < pt & qValues < qt) != 0) {
                        message('No SV while p-values and q-values have results')
                        layout(matrix(rep(
                                c(1, 1, 2, 2, 3, 3, 4, 4, 5), 9
                        ), 9, 9, byrow = TRUE))
                        par(mar = c(3, 5, 1, 1))
                        data <- data[pValues < pt & qValues < qt, ]
                        signal <-
                                signal[pValues < pt &
                                               qValues < qt, ]
                        error <-
                                error[pValues < pt & qValues < qt, ]
                        zlim <- range(c(data, signal, error))
                        
                        image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        image(
                                t(signal),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        signal
                                ) - 1)),
                                labels = colnames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        signal
                                ) - 1)),
                                labels = rownames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        image(
                                t(error),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error
                                ) - 1)),
                                labels = colnames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        error
                                ) - 1)),
                                labels = rownames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                        li <-
                                list(data, pValues < pt &
                                             qValues < qt)
                        names(li) <- c('data', 'pqvalues')
                        return(li)
                }
                else{
                        message('No SV while p-values and q-values have no results')
                        layout(matrix(rep(
                                c(1, 1, 1, 2, 2, 3, 3, 4), 8
                        ), 8, 8, byrow = TRUE))
                        par(mar = c(3, 6, 2, 3))
                        zlim <- range(c(data, signal, error))
                        
                        image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 3, 2, 1))
                        image(
                                t(signal),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        signal
                                ) - 1)),
                                labels = colnames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 3, 2, 1))
                        image(
                                t(error),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error
                                ) - 1)),
                                labels = colnames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        par(mar = c(3, 0, 2, 3))
                        plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                }
        } else{
                if (pqvalues == "anova" & sum(pValues < pt & qValues < qt) != 0) {
                        message('Have SVs while p-values and q-values have results')
                        layout(matrix(rep(
                                c(1, 1, 2, 2, 3, 3, 4, 4, 5), 9
                        ), 9, 9, byrow = TRUE))
                        par(mar = c(3, 5, 1, 1))
                        data <- data[pValues < pt & qValues < qt, ]
                        signal <-
                                signal2[pValues < pt &
                                                qValues < qt, ]
                        error <-
                                error2[pValues < pt &
                                               qValues < qt, ]
                        zlim <- range(c(data, signal, error))
                        
                        image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        image(
                                t(signal),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        signal
                                ) - 1)),
                                labels = colnames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        signal
                                ) - 1)),
                                labels = rownames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        image(
                                t(error),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error
                                ) - 1)),
                                labels = colnames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        error
                                ) - 1)),
                                labels = rownames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                        li <-
                                list(data, pValues < pt &
                                             qValues < qt)
                        names(li) <- c('data', 'pqvalues')
                        return(li)
                }
                else if (pqvalues == "anova") {
                        message('Have SVs while p-values and q-values have no results')
                        layout(matrix(rep(
                                c(1, 1, 1, 2, 2, 3, 3, 4), 8
                        ), 8, 8, byrow = TRUE))
                        par(mar = c(3, 6, 2, 3))
                        zlim <- range(c(data, signal2, error2))
                        
                        image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 3, 2, 1))
                        image(
                                t(signal2),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (
                                        ncol(signal2) - 1
                                )),
                                labels = colnames(signal2),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 3, 2, 1))
                        image(
                                t(error2),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error2
                                ) - 1)),
                                labels = colnames(error2),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        par(mar = c(3, 0, 2, 3))
                        plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                }
                else if (pqvalues == "sv" &
                         sum(pValuesSv < pt &
                             qValuesSv < qt) != 0) {
                        message('SVs corrected while p-values and q-values have results')
                        layout(matrix(rep(
                                c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6), 11
                        ), 11, 11, byrow = TRUE))
                        par(mar = c(3, 4, 2, 1))
                        data <-
                                data[pValuesSv < pt &
                                             qValuesSv < qt, ]
                        signal <- signal[pValuesSv < pt &
                                                 qValuesSv < qt, ]
                        batch <-
                                batch[pValuesSv < pt &
                                              qValuesSv < qt, ]
                        error <-
                                error[pValuesSv < pt &
                                              qValuesSv < qt, ]
                        datacor <-
                                datacor[pValuesSv < pt &
                                                qValuesSv < qt, ]
                        zlim <-
                                range(c(data, signal, batch, error, datacor))
                        
                        image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 1
                        )
                        
                        image(
                                t(signal),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        signal
                                ) - 1)),
                                labels = colnames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        signal
                                ) - 1)),
                                labels = rownames(signal),
                                cex.axis = 0.618,
                                las = 1
                        )
                        
                        image(
                                t(batch),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-batch',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        batch
                                ) - 1)),
                                labels = colnames(batch),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        batch
                                ) - 1)),
                                labels = rownames(batch),
                                cex.axis = 0.618,
                                las = 1
                        )
                        
                        image(
                                t(error),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error
                                ) - 1)),
                                labels = colnames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        error
                                ) - 1)),
                                labels = rownames(error),
                                cex.axis = 0.618,
                                las = 1
                        )
                        
                        image(
                                t(datacor),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-corrected',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (
                                        ncol(datacor) - 1
                                )),
                                labels = colnames(datacor),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (
                                        nrow(datacor) - 1
                                )),
                                labels = rownames(datacor),
                                cex.axis = 0.618,
                                las = 1
                        )
                        
                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        par(mar = c(3, 0, 2, 3))
                        plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                        li <-
                                list(datacor,
                                     data,
                                     pValuesSv < pt &
                                             qValuesSv < qt)
                        names(li) <-
                                c('datacor', 'data', 'pqvalues')
                        return(li)
                }
                else{
                        message('SVs corrected while p-values and q-values have no results')
                        layout(matrix(rep(
                                c(1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6), 13
                        ), 13, 13, byrow = TRUE))
                        par(mar = c(3, 6, 2, 3))
                        zlim <-
                                range(c(signal, data, batch, error, datacor))
                        
                        image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 3, 2, 1))
                        image(
                                t(signal),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        signal
                                ) - 1)),
                                labels = colnames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 3, 2, 1))
                        image(
                                t(batch),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-batch',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        batch
                                ) - 1)),
                                labels = colnames(batch),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 3, 2, 1))
                        image(
                                t(error),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error
                                ) - 1)),
                                labels = colnames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        par(mar = c(3, 6, 2, 3))
                        image(
                                t(datacor),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-corrected',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        axis(
                                1,
                                at = seq(0, 1, 1 / (
                                        ncol(datacor) - 1
                                )),
                                labels = colnames(datacor),
                                cex.axis = 0.618,
                                las = 2
                        )
                        
                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        par(mar = c(3, 0, 2, 3))
                        plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                }
        }
}

svaanno <-
        function(raw,
                 lv,
                 polarity = "positive",
                 projectname = "test") {
                table <- MAITbuilder(
                        data = raw$dataCorrected,
                        spectraID = NULL,
                        masses = raw$mz,
                        rt = raw$rt,
                        spectraEstimation = TRUE,
                        classes = lv,
                        rtRange = 0.2,
                        corThresh = 0.7
                )
                importMAIT <- Biotransformations(
                        MAIT.object = table,
                        adductAnnotation = TRUE,
                        peakPrecision = 0.005,
                        adductTable = NULL
                )
                importMAIT <- identifyMetabolites(
                        MAIT.object = importMAIT,
                        peakTolerance = 0.005,
                        polarity = polarity,
                        projectname = projectname
                )
                return(importMAIT)
        }