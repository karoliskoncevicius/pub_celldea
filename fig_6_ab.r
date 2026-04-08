#--- init ----------------------------------------------------------------------

library(annmatrix)
library(basetheme)
library(kklibrary)
library(EpiDISH)

X <- readRDS("data/johannson/3_annmatrix.rds")


#--- get cell counts -----------------------------------------------------------

X <- X[,X$age > 20]

C <- epidish(X, centDHSbloodDMC.m)$estF
prop <- rowSums(C[,c("CD4T", "CD8T", "NK")])


#--- plot ----------------------------------------------------------------------

cg <- "cg18263572"

cols <- c("mediumvioletred", "darkgray", "dodgerblue3")[cut(prop, c(0, 0.2, 0.32, 1))]

pdf("fig_6_ab.pdf", width = 14/2.54, height = 4/2.54, pointsize = 8)

par(mfrow = c(1,4), mar = c(3,3,3,0), oma = c(0,0,0,1), tck = -0.02, mgp = c(2,0.5,0), las = 1)

vals <- X[cg,] * 100
plot(X$age, vals, col = cols, cex = 0.5, xlab = "Age (years)", ylab = "", xaxt = "n", yaxt = "n", xlim = c(20,100), ylim = c(20,60))
axis(1, at = seq(20,100,20))
axis(2, at = c(20,30,40,50,60))
for(cl in unique(cols)) {
  abline(lm(vals ~ X$age, subset = cols == cl), col = cl)
}
title(cg, line = 0.5)
title(ylab = "Modification (%)", line = 1.7)
mtext("a", 3, adj = -0.2, line = 1.7, font = 2)

vals <- residuals(lm(X[cg,] ~ X$age))
plot(X$age, vals, col = cols, cex = 0.5, xlab = "Age (years)", ylab = "", ylim = c(-0.3, 0.3), xaxt = "n", yaxt = "n", xlim = c(20,100))
axis(1, at = seq(20,100,20))
axis(2, at = c(-0.2,0,0.2))
for(cl in unique(cols)) {
  abline(lm(vals ~ X$age, subset = cols == cl), col = cl)
}
title(ylab = "Residual", line = 1.7)
mtext("b", 3, adj = -0.2, line = 1.7, font = 2)

vals <- residuals(lm(X[cg,] ~ X$age + C))
plot(X$age, vals, col = cols, cex = 0.5, xlab = "Age (years)", ylab = "", ylim = c(-0.3, 0.3), xaxt = "n", yaxt = "n", xlim = c(20,100))
axis(1, at = seq(20,100,20))
axis(2, at = c(-0.2,0,0.2))
for(cl in unique(cols)) {
  abline(lm(vals ~ X$age, subset = cols == cl), col = cl)
}
title(ylab = "Residual", line = 1.7)

vals <- residuals(lm(X[cg,] ~ X$age * C))
plot(X$age, vals, col = cols, cex = 0.5, xlab = "Age (years)", ylab = "", ylim = c(-0.3, 0.3), xaxt = "n", yaxt = "n", xlim = c(20,100))
axis(1, at = seq(20,100,20))
axis(2, at = c(-0.2,0,0.2))
for(cl in unique(cols)) {
  abline(lm(vals ~ X$age, subset = cols == cl), col = cl)
}
title(ylab = "Residual", line = 1.7)

legend("bottom", legend = c("<20 %", "20 % - 30 %", ">30 %"), text.width = strwidth(c("<20 %", "20 % - 30 %", ">30 %")), col = c("mediumvioletred", "darkgray", "dodgerblue3"),
       pch = 19, bty = 'n', inset = c(0, 1), xpd = NA, horiz = TRUE, title = "T lymphocyte %")

invisible(dev.off())
