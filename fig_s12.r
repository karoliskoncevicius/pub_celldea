#--- init ----------------------------------------------------------------------

library(annmatrix)
library(basetheme)
library(kklibrary)
library(EpiDISH)
library(DEAfree)

X <- readRDS("data/salas/3_annmatrix.rds")


#--- get cell proportions ------------------------------------------------------

X <- X[,X$celltype == "mix"]
C <- epidish(X, centDEAfree12CT.m)$estF * 100


#--- plot ----------------------------------------------------------------------

real <- cbind(neu = X$perc_neu, mono = X$perc_mono, eos = X$perc_eos,
              nk = X$perc_nk, bn = X$perc_bn, bm = X$perc_bm, treg = X$perc_treg,
              tcd4n = X$perc_tcd4n, tcd4m = X$perc_tcd4m, tcd8 = X$perc_tcd8) * 100

C <- cbind(C, tcd8 = C[,"tcd8n"] + C[,"tcd8m"])
C <- C[,colnames(real)]

N <- c(neu = "Neutrophils", mono = "Monocytes", eos = "Eosinophils", nk = "Natural killers", bn = "B naive", bm = "B memory", treg = "T regulatory", tcd4n = "T CD4 naive", tcd4m = "T CD4 memory", tcd8 = "T CD8")


pdf("fig_s12.pdf", width = 16/2.54, height = 7/2.54, pointsize = 8)

par(mfrow = c(2,5), mar = c(3,3,3,0), oma = c(0,0,0,1), tck = -0.02, mgp = c(2,0.5,0), las = 1)

for(i in 1:ncol(real)) {
  plot.new()
  plot.window(xlim = c(0,30), ylim = c(0,30))
  abline(0,1, col = "gray", lty = 2)
  points(real[,i], C[,i], col = cell2col(colnames(C)[i]))
  axis(1)
  axis(2)
  title(N[colnames(C)[i]], line = 0.5)
  title(xlab = "Real proportion (%)")
  title(ylab = "Estimated proportion (%)", line = 1.65)
  box()
  rmse <- round(sqrt(mean((C[,i] - real[,i])^2)),2)
  cor  <- round(cor(C[,i], real[,i]), 2)
  legend("topleft", legend = paste0(c("RMSE = ", "r = "), c(rmse, cor)), inset = c(-0.075,0), bty = "n")
}

invisible(dev.off())
