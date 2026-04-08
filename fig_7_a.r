#--- init ----------------------------------------------------------------------

library(annmatrix)
library(kklibrary)
library(matrixTests)
library(EpiDISH)
library(DEAfree)
library(methylCIPHER)
suppressPackageStartupMessages(library(inops))

X <- readRDS("data/hannon_scz/1_normalized.rds")

cent1 <- centDHSbloodDMC.m
cent2 <- cent12CT450k.m
cent3 <- centDEAfree12CT.m

# rename and reorder
newnames <- c(Neutro = "Neu", Mono = "Mono", Eosino = "Eos", NK = "NK", B = "Bcell", CD4T = "TCD4", CD8T = "TCD8")
colnames(cent1) <- newnames[colnames(cent1)]
cent1           <- cent1[,newnames]

newnames <- c(Neu = "Neu", Mono = "Mono", Eos = "Eos", Baso = "Baso", NK = "NK", Bnv = "Bn", Bmem = "Bm", Treg = "Treg", CD4Tnv = "TCD4n", CD4Tmem = "TCD4m", CD8Tnv = "TCD8n", CD8Tmem = "TCD8m")
colnames(cent2) <- newnames[colnames(cent2)]
cent2           <- cent2[,newnames]

cent3 <- cent3[,tolower(newnames)]
colnames(cent3) <- newnames


#--- clean ---------------------------------------------------------------------

X <- X[,!is.na(X$sex) & !is.na(X$age) & X$age < 100 & X$cohort %in% c("gse147221", "gse80417", "gse84727")]

# smoking score
X$smokingscore <- calcSmokingMcCartney(t(X))


#--- celltype proportion differences -------------------------------------------

C1 <- epidish(X, cent1)$estF
C2 <- epidish(X, cent2)$estF
C3 <- epidish(X, cent3)$estF

plotCells <- function(X, C) {

  res <- list()
  for(cell in colnames(C)) {
    res[[cell]] <- lm(C[,cell] ~ X$diagnosis + X$age + X$sex + X$cohort + X$smokingscore)
  }
  beta <- mapply(function(x) summary(x)$coefficients[2,1], res) * 100
  pval <- mapply(function(x) summary(x)$coefficients[2,4], res)

  plot.new()
  plot.window(xlim = c(-5,15), ylim = c(length(beta)+0.5, -0.5), xaxs = "i", yaxs = "i")

  text(-5, 1:length(beta), names(beta), xpd = NA, font = 2, pos = 2)

  abline(v = 0, lty = 3, col = "gray")
  points(beta, 1:length(beta), col = cell2col(tolower(names(beta))), pch = 18, cex = 3)

  axis(1, at = c(-5,0,5))
  mtext(bquote("SCZ" ~ beta ~ "coefficient"), 1, at = 0, font = 2, line = 3, cex = 0.7)

  text(8.5, 0, bquote(beta), font = 2, xpd = TRUE)
  text(12.5, 0, bquote("p-value"), font = 2, xpd = TRUE)
  text(8.5, 1:length(beta), sprintf("%0.1f", beta))
  text(12.5, 1:length(pval), formatC(pval, digits = 2))
  text(13.75, 1:length(res), ifelse(p.adjust(pval, "bonf") < 0.05, "*", ""), pos = 4, xpd = NA)

  box()

}


#--- plot ----------------------------------------------------------------------

pdf("fig_7_a.pdf", width = 14/2.54, height = 5/2.54, pointsize = 8)

par(mfrow = c(1,3), mar = c(4,3,3,1))

plotCells(X, C1)
title("EpiDISH-7", line = 0.5)
mtext("a", 3, line = 1.9, font = 2, adj = -0.13)
plotCells(X, C2)
title("EpiDISH-12", line = 0.5)
plotCells(X, C3)
title("DEA-free", line = 0.5)

invisible(dev.off())
