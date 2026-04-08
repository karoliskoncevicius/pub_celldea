#--- init ----------------------------------------------------------------------

library(annmatrix)
library(basetheme)
library(kklibrary)
library(EpiDISH)
library(DEAfree)

X <- readRDS("data/rhead/3_annmatrix.rds")
C <- readRDS("data/rhead/cellcountsepidish.rds")

C <- C[colnames(X),] * 100
stopifnot(all.equal(rownames(C), colnames(X)))


#--- get cell proportions ------------------------------------------------------

C2 <- epidish(X, centDEAfree12CT.m)$estF * 100


#--- plot ----------------------------------------------------------------------

plot_purity <- function(age, pur, col) {
  plot(age, pur, col = col, ylim = c(60,100), xlab = "Age (years)", ylab = "")
  title(ylab = "Estimated purity (%)", line = 1.6)

  abline(lm(pur ~ age), col = "red2", lty = 2)

  txt <- paste0(c("r = ", "p = "), c(round(cor(pur, age), 2), formatC(cor.test(pur, age)$p.value, digits = 2)))
  legend("bottomleft", legend = txt, text.col = "red2", bty = "n", inset = c(-0.075,0))
}


pdf("fig_5_de.pdf", width = 14/2.54, height = 4/2.54, pointsize = 8)

par(mfrow = c(1,4), mar = c(3,3,3,0), oma = c(0,0,0,1), tck = -0.02, mgp = c(2,0.5,0), las = 1)

# tcd4n epidish
plot_purity(X$age[X$celltype == "tcd4n"], C[X$celltype == "tcd4n", "tcd4n"], cell2col("tcd4n"))
title("EpiDISH-12", adj = 0, line = 0.5)
mtext("d", 3, adj = -0.2, line = 1.7, font = 2)

# tcd4m epidish
plot_purity(X$age[X$celltype == "tcd4m"], C[X$celltype == "tcd4m", "tcd4m"], cell2col("tcd4m"))


# tcd4n deafree
plot_purity(X$age[X$celltype == "tcd4n"], C2[X$celltype == "tcd4n", "tcd4n"], cell2col("tcd4n"))
title("DEA-free", adj = 0, line = 0.5)
mtext("e", 3, adj = -0.2, line = 1.7, font = 2)

# tcd4m deafree
plot_purity(X$age[X$celltype == "tcd4m"], C2[X$celltype == "tcd4m", "tcd4m"], cell2col("tcd4m"))

legend("bottomright", legend = c("TCD4n", "TCD4m"), col = cell2col(c("tcd4n", "tcd4m")), pch = 19, bty = 'n', inset = c(0, 1.025), xpd = NA, horiz = TRUE)


invisible(dev.off())
