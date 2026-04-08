#--- init ----------------------------------------------------------------------

library(annmatrix)
library(basetheme)
library(kklibrary)
library(EpiDISH)
suppressPackageStartupMessages(library(FlowSorted.BloodExtended.EPIC))

X     <- readRDS("data/rhead/3_annmatrix.rds")
dfage <- readRDS("data/diffage_cgs.rds")

X <- X[rownames(dfage),]
stopifnot(all.equal(rownames(X), rownames(dfage)))

# loading the manually-saved panels for minfi
minfi6  <- read.csv("data/Minfi/minfi_6ct.csv", row.names = 1)
minfi12 <- read.csv("data/Minfi/minfi_12ct.csv", row.names = 1)


#--- plot ----------------------------------------------------------------------

# so that celltypes do not overlap each other on the plot
X$celltype <- droplevels(factor(X$celltype, levels = c("neu", "mono", "bcell", "bn", "bm", "nk", "tcell", "tcd4n", "tcd4m", "tcd8n", "tcd8m")))
X <- X[, sample(1:ncol(X))]

# NOTE: manually selected cytosines from epidish reference
cgs <- c("cg07728865", "cg15163862")

pdf("fig_5_abc.pdf", width = 14/2.54, height = 4/2.54, pointsize = 8)

par(mfrow = c(1,4), mar = c(3,3,3,0), oma = c(0,0,0,1), tck = -0.02, mgp = c(2,0.5,0), las = 1)

# overlaps
perc <- pval <- numeric()

perc[1] <- mean((dfage$category == "cell-dea")[rownames(dfage) %in% rownames(minfi6)], na.rm = TRUE) * 100
perc[2] <- mean((dfage$category == "cell-dea")[rownames(dfage) %in% rownames(minfi12)], na.rm = TRUE) * 100
perc[3] <- mean((dfage$category == "cell-dea")[rownames(dfage) %in% IDOLOptimizedCpGs450klegacy], na.rm = TRUE) * 100
perc[4] <- mean((dfage$category == "cell-dea")[rownames(dfage) %in% IDOLOptimizedCpGsBloodExtended450k], na.rm = TRUE) * 100
perc[5] <- mean((dfage$category == "cell-dea")[rownames(dfage) %in% rownames(centDHSbloodDMC.m)], na.rm = TRUE) * 100
perc[6] <- mean((dfage$category == "cell-dea")[rownames(dfage) %in% rownames(cent12CT450k.m)], na.rm = TRUE) * 100
perc[7] <- mean((dfage$category == "cell-dea")[rownames(dfage) %in% rownames(centUniLIFE.m)], na.rm = TRUE) * 100

pval[1] <- fisher.test(table(dfage$category == "cell-dea", rownames(dfage) %in% rownames(minfi6)))$p.value
pval[2] <- fisher.test(table(dfage$category == "cell-dea", rownames(dfage) %in% rownames(minfi12)))$p.value
pval[3] <- fisher.test(table(dfage$category == "cell-dea", rownames(dfage) %in% IDOLOptimizedCpGs450klegacy))$p.value
pval[4] <- fisher.test(table(dfage$category == "cell-dea", rownames(dfage) %in% IDOLOptimizedCpGsBloodExtended450k))$p.value
pval[5] <- fisher.test(table(dfage$category == "cell-dea", rownames(dfage) %in% rownames(centDHSbloodDMC.m)))$p.value
pval[6] <- fisher.test(table(dfage$category == "cell-dea", rownames(dfage) %in% rownames(cent12CT450k.m)))$p.value
pval[7] <- fisher.test(table(dfage$category == "cell-dea", rownames(dfage) %in% rownames(centUniLIFE.m)))$p.value

names(perc) <- c("Minfi-6", "Minfi-12", "IDOL-6", "IDOL-12", "EpiDISH-7", "EpiDISH-12", "EpiDISH-19")


par(mar = c(3,5.5,3,0))
xs <- barplot(perc, xlim = c(0,100), col = "#E8E8E8", horiz = TRUE, las = 1)
title(xlab = "Cell-DEA cytosines (%)", line = 2)
text(perc, xs, paste("p = ", formatC(pval, digit = 2)), pos = 4, col = "red2")

mtext("a", 3, adj = -0.5, line = 1.7, font = 2)


# by celltype
perc <- tapply(rownames(cent12CT450k.m) %in% rownames(dfage)[dfage$category == "cell-dea"], rep(colnames(cent12CT450k.m), each = 50), mean) * 100

newnames    <- rev(c(Bnv = "Bn", Bmem = "Bm", CD4Tnv = "TCD4n", CD4Tmem = "TCD4m", CD8Tnv = "TCD8n", CD8Tmem = "TCD8m", Treg = "Treg", NK = "NK", Neu = "Neu", Mono = "Mono", Eos = "Eos", Baso = "Baso"))
names(perc) <- newnames[names(perc)]  # rename
perc        <- perc[newnames]         # reorder

par(mar = c(3,3,3,0))
barplot(perc, xlim = c(0,100), col = cell2col(tolower(names(perc))), horiz = TRUE, las = 1)
title(xlab = "Cell-DEA cytosines (%)", line = 2)

mtext("b", 3, adj = -0.25, line = 1.7, font = 2)


# example cytosines
par(mar = c(3,3,3,0))

for(cg in cgs) {
  plot(X$age, X[cg,] * 100, col = cell2col(X$celltype), pch = 19, cex = 0.8, ylim = c(0,100), xlab = "", ylab = "", yaxt = "n", yaxs = "i")
  abline(h = c(cent12CT450k.m[cg,c("CD4Tnv", "CD4Tmem", "Mono")], sum(cent12CT450k.m[cg,c("Bnv","Bmem")])) * 100, col = cell2col(c("tcd4n", "tcd4m", "mono", "bcell")))
  title(cg, line = 0.5)
  title(ylab = "Modification (%)", line = 1.6)
  axis(2, las = 2)
  title(xlab = "Age (years)")
  if(cg == cgs[1]) {
     mtext("c", 3, adj = -0.2, line = 1.7, font = 2)
  }
}

legend("bottom", legend = c("Mono", "TCD4n", "TCD4m"), text.width = strwidth(c("Mono", "TCD4n", "TCD4m")), col = cell2col(c("mono", "tcd4n", "tcd4m")),
       pch = 19, bty = 'n', inset = c(0, 1.1), xpd = NA, horiz = TRUE)


invisible(dev.off())
