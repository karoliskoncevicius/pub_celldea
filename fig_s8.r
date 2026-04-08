#--- init ----------------------------------------------------------------------

library(annmatrix)
library(basetheme)
library(kklibrary)
library(EpiDISH)
library(DEAfree)

X  <- readRDS("data/blueprint/3_annmatrix.rds")
C1 <- readRDS("data/blueprint/cellcountsepidish.rds")
C2 <- readRDS("data/blueprint/cellcountsidol.rds")
C3 <- readRDS("data/blueprint/cellcountsminfi.rds")

C1 <- C1[colnames(X),] * 100
C2 <- C2[colnames(X),] * 100
C3 <- C3[colnames(X),] * 100

stopifnot(all.equal(rownames(C1), colnames(X)))
stopifnot(all.equal(rownames(C2), colnames(X)))
stopifnot(all.equal(rownames(C3), colnames(X)))


#--- get cell proportions ------------------------------------------------------

C4 <- epidish(X, centDEAfree12CT.m)$estF * 100


#--- plot ----------------------------------------------------------------------

plot_purity <- function(age, pur, col, side = "bottomleft") {
  plot(age, pur, col = col, ylim = c(50,100), xlab = "Age (years)", ylab = "")
  title(ylab = "Estimated purity (%)", line = 1.6)

  abline(lm(pur ~ age), col = "red2", lty = 2)

  txt <- paste0(c("r = ", "p = "), c(round(cor(pur, age), 2), formatC(cor.test(pur, age)$p.value, digits = 2)))
  legend(side, legend = txt, text.col = "red2", bty = "n", inset = c(-0.075,0))
}


pdf("fig_s8.pdf", width = 14/2.54, height = 4/2.54, pointsize = 8)

par(mfrow = c(1,4), mar = c(3,3,3,0), oma = c(0,0,0,1), tck = -0.02, mgp = c(2,0.5,0), las = 1)

# tcd4n epidish
plot_purity(X$age[X$celltype == "tcd4n"], C1[X$celltype == "tcd4n", "tcd4n"], cell2col("tcd4n"))
title("EpiDISH-12", adj = 0, line = 0.5)

# tcd4n idol
plot_purity(X$age[X$celltype == "tcd4n"], C2[X$celltype == "tcd4n", "tcd4n"], cell2col("tcd4n"))
title("IDOL-12", adj = 0, line = 0.5)

# tcd4n minfi
plot_purity(X$age[X$celltype == "tcd4n"], C3[X$celltype == "tcd4n", "tcd4n"], cell2col("tcd4n"), "topleft")
title("Minfi-12", adj = 0, line = 0.5)

# tcd4n deafree
plot_purity(X$age[X$celltype == "tcd4n"], C4[X$celltype == "tcd4n", "tcd4n"], cell2col("tcd4n"))
title("DEA-free", adj = 0, line = 0.5)


invisible(dev.off())
