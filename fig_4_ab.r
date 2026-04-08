#--- init ----------------------------------------------------------------------

library(readxl)
library(annmatrix)
library(kklibrary)
suppressPackageStartupMessages(library(inops))

dfage <- readRDS("data/diffage_cgs.rds")
ewas  <- as.data.frame(read_excel("data/EWAS/elife-58430-supp1-v2.xlsx", sheet = 3, skip = 3))


#--- clean ---------------------------------------------------------------------

rownames(ewas) <- ewas[,1]
ewas <- ewas[,c(7,9)]
colnames(ewas) <- c("diff", "pvalue")

cgs   <- intersect(rownames(ewas), rownames(dfage))
ewas  <- ewas[cgs,]
dfage <- dfage[cgs,]


#--- plot ----------------------------------------------------------------------

myplot <- function(age, diff, col, ...) {
  plot(age, diff, col = adjustcolor(col, 0.3), pch = 19, cex = 0.6, xlab = "", ylab = "", axes = FALSE, ...)
  abline(h = 0, col = "gray", lty = 2)
  abline(v = 0, col = "gray", lty = 2)
  axis(2, at = c(-2,0,2), las = 2)
  title(xlab = bquote(Delta ~ "Modification (%/year)"))
  title(ylab = bquote("SCZ" ~ beta ~ "coefficient"), line = 1.25)
  box()

  fish1 <- fisher.test(table(age < 0, diff > 0))
  fish2 <- fisher.test(table(age > 0, diff > 0))
  fish  <- if(fish1$estimate > fish2$estimate) fish1 else fish2
  fish$estimate <- ifelse(fish$estimate > 100, round(fish$estimate), ifelse(fish$estimate > 10, round(fish$estimate, 1), round(fish$estimate, 2)))
  fish$p.value <- ifelse(fish$p.value < 0.001, formatC(fish$p.value, format = "e", digits = 1), signif(fish$p.value, 2))
  legend("topleft", paste(c("OR =", "p ="), c(fish$estimate, fish$p.value)), text.col = "red2", bty = "n", inset = c(-0.1,-0.03))
}

pdf("fig_4_ab.pdf", width = 14/2.54, height = 4/2.54, pointsize = 8)

par(mfrow = c(1,4), mar = c(3,3,3,0), oma = c(0,0,0,1), tck = -0.02, mgp = c(2,0.5,0), las = 1)

cgs <- rownames(dfage[dfage$category == "cell-dea",])
myplot(dfage[cgs,]$ra_agingrate_tcd4n, ewas[cgs,]$diff, cell2col("tcd4n"), xlim = c(-0.6,0.6), ylim = c(-2.2,2.8))
axis(1, at = c(-0.4, 0, 0.4))
title("Cell-DEA", adj = 0, line = 0.5)
mtext("a", 3, adj = -0.2, line = 1.7, font = 2)

myplot(dfage[cgs,]$ra_agingrate_tcd4m, ewas[cgs,]$diff, cell2col("tcd4m"), xlim = c(-0.6,0.6), ylim = c(-2.2,2.8))
axis(1, at = c(-0.4, 0, 0.4))

cgs <- rownames(dfage[dfage$category != "cell-dea",])
myplot(dfage[cgs,]$ra_agingrate_tcd4n, ewas[cgs,]$diff, cell2col("tcd4n"), xlim = c(-0.6,0.6), ylim = c(-2.2,2.8))
axis(1, at = c(-0.4, 0, 0.4))
title("Non-cell-DEA", adj = 0, line = 0.5)
mtext("b", 3, adj = -0.2, line = 1.7, font = 2)

myplot(dfage[cgs,]$ra_agingrate_tcd4m, ewas[cgs,]$diff, cell2col("tcd4m"), xlim = c(-0.6,0.6), ylim = c(-2.2,2.8))
axis(1, at = c(-0.4, 0, 0.4))

legend("bottomright", legend = c("TCD4n", "TCD4m"), col = cell2col(c("tcd4n", "tcd4m")), pch = 19, bty = 'n', inset = c(0, 1.025), xpd = NA, horiz = TRUE)

invisible(dev.off())
