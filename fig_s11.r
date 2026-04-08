#--- init ----------------------------------------------------------------------

library(annmatrix)
library(basetheme)
library(kklibrary)
library(EpiDISH)

X <- readRDS("data/koncevicius/3_annmatrix.rds")
C <- readRDS("data/koncevicius/cellcountsepidish.rds")


#--- get cell proportions ------------------------------------------------------

C <- C[colnames(X),] * 100


#--- plot ----------------------------------------------------------------------

plot_purity <- function(age, pur, col, ylim = c(0,100)) {
  plot(age, pur, col = col, xlab = "Age (years)", ylab = "", xlim = c(20,70), ylim = ylim)
  abline(lm(pur ~ age), col = "red2", lty = 2)
  title(ylab = "Proportion estimate (%)", line = 2)
  cor <- cor.test(age, pur)
  txt <- paste0(c("r = ", "p = "), c(round(cor$estimate, 2), formatC(cor$p.value, digits = 2)))
  place <- ifelse(ylim[1] < 50, "topleft", "bottomleft")
  legend(place, legend = txt, text.col = "red2", bty = "n", inset = c(-0.075,0))
}

plot_cor <- function(pur1, pur2, col, ylim = c(0,100)) {
  plot(pur1, pur2, col = col, xlab = "Cytometry estimate (%)", ylab = "", xlim = ylim, ylim = ylim)
  title(ylab = "Deconvolution estimate (%)", line = 2)
  abline(0, 1, lty = 2, col = "gray")
  rmse <- sqrt(mean((pur2-pur1)^2, na.rm = TRUE))
  cor  <- cor.test(pur1, pur2)
  txt  <- paste0(c("RMSE = ", "r = "), c(round(rmse, 2), round(cor$estimate, 2)))
  legend("topleft", legend = txt, text.col = "red2", bty = "n", inset = c(-0.075,0))
}

pdf("fig_s11.pdf", width = 10/2.54, height = 13/2.54, pointsize = 8)

par(mfrow = c(4,3), mar = c(3,4,2,0), oma = c(0,0,0,1), tck = -0.02, mgp = c(2,0.5,0), las = 1)

# neutrophils
inds <- X$celltype == "neu" & !is.na(X$purity)
plot_purity(X$age[inds], X$purity[inds], cell2col("neu"), ylim = c(90,100))
title("Cytometry", line = 0.5)

plot_purity(X$age[inds], C[inds,"neu"], cell2col("neu"), ylim = c(90,100))
title("EpiDISH-12", line = 0.5)

plot_cor(X$purity[inds], C[inds,"neu"], cell2col("neu"), ylim = c(90,100))
title("Correlation", line = 0.5)


# T CD4 nv cells
inds <- X$celltype == "tcell" & !is.na(X$perc_tcd4n)
plot_purity(X$age[inds], X$perc_tcd4n[inds], cell2col("tcd4n"))
plot_purity(X$age[inds], C[inds,"tcd4n"], cell2col("tcd4n"))
plot_cor(X$perc_tcd4n[inds], C[inds,"tcd4n"], cell2col("tcd4n"))


# T CD4 mem cells
inds <- X$celltype == "tcell" & !is.na(X$perc_tcd4m)
plot_purity(X$age[inds], X$perc_tcd4m[inds], cell2col("tcd4m"))
plot_purity(X$age[inds], C[inds,"tcd4m"], cell2col("tcd4m"))
plot_cor(X$perc_tcd4m[inds], C[inds,"tcd4m"], cell2col("tcd4m"))


# T CD8 cells
inds <- X$celltype == "tcell" & !is.na(X$perc_tcd8)
plot_purity(X$age[inds], X$perc_tcd8[inds], cell2col("tcd8"))
plot_purity(X$age[inds], rowSums(C[inds,c("tcd8n", "tcd8m")]), cell2col("tcd8"))
plot_cor(X$perc_tcd8[inds], rowSums(C[inds,c("tcd8n", "tcd8m")]), cell2col("tcd8"))


invisible(dev.off())
