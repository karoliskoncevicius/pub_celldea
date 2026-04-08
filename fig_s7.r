#--- init ----------------------------------------------------------------------

library(annmatrix)
library(basetheme)
library(kklibrary)
library(EpiDISH)

X  <- readRDS("data/rhead/3_annmatrix.rds")
C1 <- readRDS("data/rhead/cellcountsidol.rds")
C2 <- readRDS("data/rhead/cellcountsminfi.rds")

C1 <- C1[colnames(X),] * 100
C2 <- C2[colnames(X),] * 100

stopifnot(all.equal(rownames(C1), colnames(X)))
stopifnot(all.equal(rownames(C2), colnames(X)))


#--- plot ----------------------------------------------------------------------

plot_purity <- function(age, pur, col, side = "bottomleft") {
  plot(age, pur, col = col, ylim = c(50,100), xlab = "Age (years)", ylab = "")
  title(ylab = "Estimated purity (%)", line = 1.6)

  abline(lm(pur ~ age), col = "red2", lty = 2)

  txt <- paste0(c("r = ", "p = "), c(round(cor(pur, age), 2), formatC(cor.test(pur, age)$p.value, digits = 2)))
  legend(side, legend = txt, text.col = "red2", bty = "n", inset = c(-0.075,0))
}


pdf("fig_s7.pdf", width = 14/2.54, height = 4/2.54, pointsize = 8)

par(mfrow = c(1,4), mar = c(3,3,3,0), oma = c(0,0,0,1), tck = -0.02, mgp = c(2,0.5,0), las = 1)

# idol
plot_purity(X$age[X$celltype == "tcd4n"], C1[X$celltype == "tcd4n", "tcd4n"], cell2col("tcd4n"))
title("IDOL-12", adj = 0, line = 0.5)
plot_purity(X$age[X$celltype == "tcd4m"], C1[X$celltype == "tcd4m", "tcd4m"], cell2col("tcd4m"))

# minfi
plot_purity(X$age[X$celltype == "tcd4n"], C2[X$celltype == "tcd4n", "tcd4n"], cell2col("tcd4n"), "topleft")
title("Minfi-12", adj = 0, line = 0.5)
plot_purity(X$age[X$celltype == "tcd4m"], C2[X$celltype == "tcd4m", "tcd4m"], cell2col("tcd4m"), "topleft")


invisible(dev.off())
