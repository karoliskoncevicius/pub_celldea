#--- init ----------------------------------------------------------------------

library(annmatrix)
library(basetheme)
library(kklibrary)
library(EpiDISH)
library(DEAfree)

X <- readRDS("data/blueprint2/3_annmatrix.rds")
X <- X[,!is.na(X$facs_purity) & X$celltype == "tcd4"]


#--- get cell proportions ------------------------------------------------------

C1 <- epidish(X, centDHSbloodDMC.m)$estF
C2 <- epidish(X, cent12CT450k.m)$estF
C3 <- epidish(X, centUniLIFE.m)$estF
C0 <- epidish(X, centDEAfree12CT.m)$estF


#--- plot ----------------------------------------------------------------------

plot_purity <- function(age, pur, col) {
  plot(age, pur, ylim = c(70,100), xlim = c(0,80), col = col, xlab = "Age (years)", yaxt = "n", ylab = "", las = 1)
  abline(lm(pur ~ age), col = "red2", lty = 2)

  txt <- paste0(c("r = ", "p = "), c(round(cor(pur, age), 2), formatC(cor.test(pur, age)$p.value, digits = 2)))
  legend("bottomleft", legend = txt, text.col = "red2", bty = "n", inset = c(-0.075,0))
}


pdf("fig_5_f.pdf", width = 14/2.54, height = 4/2.54, pointsize = 8)

par(mfrow = c(1,5), mar = c(5,0,3,1), oma = c(0,4,0,0), tck = -0.02, mgp = c(2,0.5,0))

# cytometry
plot_purity(X$age, X$facs_purity * 100, cell2col("tcd4"))
axis(2, las = 1, at = c(70,80,90,100))
title("Cytometry", line = 0.5)
mtext("f", 3, line = 1.5, font = 2, adj = -0.3)

title(ylab = "Estimated purity (%)", line = 2.5, outer = TRUE)


# decomposition estimates
plot_purity(X$age, C1[,"CD4T"] * 100, cell2col("tcd4"))
title("EpiDISH-7", line = 0.5)

plot_purity(X$age, rowSums(C2[,c("CD4Tnv", "CD4Tmem", "Treg")]) * 100, cell2col("tcd4"))
title("EpiDISH-12", line = 0.5)

plot_purity(X$age, rowSums(C3[,c("CD4T", "aCD4Tnv", "aCD4Tmem", "aTreg")]) * 100, cell2col("tcd4"))
title("EpiDISH-19", line = 0.5)

plot_purity(X$age, rowSums(C0[,c("tcd4n", "tcd4m", "treg")]) * 100, cell2col("tcd4"))
title("DEA-free", line = 0.5)


invisible(dev.off())
