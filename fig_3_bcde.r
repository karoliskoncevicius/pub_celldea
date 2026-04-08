#--- init ----------------------------------------------------------------------

library(readxl)
library(annmatrix)
library(kklibrary)
library(matrixTests)
library(EpiDISH)
library(DEAfree)
library(methylCIPHER)

X    <- readRDS("data/koncevicius/3_annmatrix.rds")
dfa  <- readRDS("data/diffage_cgs.rds")

ewasS <- as.data.frame(read_excel("data/EWAS/elife-58430-supp1-v2.xlsx", sheet = 3, skip = 3))

ewasB <- as.data.frame(read_excel("data/EWAS/41366_2018_BFijo2017269_MOESM66_ESM.xlsx", sheet = 3))
ewasB <- ewasB[!is.na(ewasB[,1]) & (startsWith(ewasB[,1], "cg") | startsWith(ewasB[,1], "ch")),]

ewasT <- as.data.frame(read_excel("data/EWAS/pone.0063812.s005.xls", skip = 2))
ewasT <- ewasT[!is.na(ewasT[,1]) & nchar(ewasT[,1]) < 15,]
ewasT[,1] <- sub("\\*", "", ewasT[,1])


#--- clean ---------------------------------------------------------------------

# get cell type proportions
C <- epidish(X, centDEAfree12CT.m)$estF

# select differences and p-values from EWAS
rownames(ewasS) <- ewasS[,1]
ewasS <- ewasS[,c(7,9)]
colnames(ewasS) <- c("diff", "pvalue")

rownames(ewasB) <- ewasB[,1]
ewasB <- ewasB[,c(12,13)]
colnames(ewasB) <- c("diff", "pvalue")

rownames(ewasT) <- ewasT[,1]
ewasT <- ewasT[,c(6,7)]
colnames(ewasT) <- c("diff", "pvalue")

# split cell types
NS <- X[rownames(X) %in% rownames(ewasS), X$celltype == "neu"]
TS <- X[rownames(X) %in% rownames(ewasS), X$celltype != "neu"]
NB <- X[rownames(X) %in% rownames(ewasB), X$celltype == "neu"]
TB <- X[rownames(X) %in% rownames(ewasB), X$celltype != "neu"]
NT <- X[rownames(X) %in% rownames(ewasT), X$celltype == "neu"]
TT <- X[rownames(X) %in% rownames(ewasT), X$celltype != "neu"]

# align ewas tables
ewasS <- ewasS[rownames(NS),]
ewasB <- ewasB[rownames(NB),]
ewasT <- ewasT[rownames(NT),]


#--- results -------------------------------------------------------------------

model  <- model.matrix(~ NS$diagnosis + NS$age + NS$sex + NS$bmi + NS$smoking)
model0 <- model.matrix(~ NS$age + NS$sex + NS$bmi + NS$smoking)
resNS  <- row_lm_f(NS, model, model0)

model  <- model.matrix(~ TS$diagnosis + TS$age + TS$sex + TS$bmi + TS$smoking + C[colnames(TS),c("tcd4n", "tcd4m", "tcd8n", "tcd8m", "treg")])
model0 <- model.matrix(~ TS$age + TS$sex + TS$bmi + TS$smoking + C[colnames(TS),c("tcd4n", "tcd4m", "tcd8n", "tcd8m", "treg")])
resTS   <- row_lm_f(TS, model, model0)

model  <- model.matrix(~ NB$bmi + NB$age + NB$sex + NB$diagnosis + NB$smoking)
model0 <- model.matrix(~ NB$age + NB$sex + NB$diagnosis + NB$smoking)
resNB  <- row_lm_f(NB, model, model0)

model  <- model.matrix(~ TB$bmi + TB$age + TB$sex + TB$diagnosis + TB$smoking + C[colnames(TB),c("tcd4n", "tcd4m", "tcd8n", "tcd8m", "treg")])
model0 <- model.matrix(~ TB$age + TB$sex + TB$diagnosis + TB$smoking + C[colnames(TB),c("tcd4n", "tcd4m", "tcd8n", "tcd8m", "treg")])
resTB  <- row_lm_f(TB, model, model0)

model  <- model.matrix(~ NT$smoking + NT$age + NT$sex + NT$diagnosis + NT$bmi)
model0 <- model.matrix(~ NT$age + NT$sex + NT$diagnosis + NT$bmi)
resNT  <- row_lm_f(NT, model, model0)

model  <- model.matrix(~ TT$smoking + TT$age + TT$sex + TT$diagnosis + TT$bmi + C[colnames(TT),c("tcd4n", "tcd4m", "tcd8n", "tcd8m", "treg")])
model0 <- model.matrix(~ TT$age + TT$sex + TT$diagnosis + TT$bmi + C[colnames(TT),c("tcd4n", "tcd4m", "tcd8n", "tcd8m", "treg")])
resTT  <- row_lm_f(TT, model, model0)


# add diff-age info
resNS$diffage <- rownames(resNS) %in% rownames(dfa)[dfa$category == "cell-dea"]
resTS$diffage <- rownames(resTS) %in% rownames(dfa)[dfa$category == "cell-dea"]
resNB$diffage <- rownames(resNB) %in% rownames(dfa)[dfa$category == "cell-dea"]
resTB$diffage <- rownames(resTB) %in% rownames(dfa)[dfa$category == "cell-dea"]
resNT$diffage <- rownames(resNT) %in% rownames(dfa)[dfa$category == "cell-dea"]
resTT$diffage <- rownames(resTT) %in% rownames(dfa)[dfa$category == "cell-dea"]


#--- plot ----------------------------------------------------------------------


pdf("fig_3_bcde.pdf", width = 14/2.54, height = 8/2.54, pointsize = 8)

par(mfrow = c(2,4), mar = c(3,3,3,0), oma = c(0,0,0,1), tck = -0.02, mgp = c(2,0.5,0), las = 1)


# schizophrenia neutrophils
fdrcut <- mean(c(max(resNS$pvalue[p.adjust(resNS$pvalue, "fdr") < 0.05]), min(resNS$pvalue[p.adjust(resNS$pvalue, "fdr") > 0.05])))

plot.new()
plot.window(xlim = range(resNS$beta.1 * 100), ylim = c(0,6))
abline(h = -log10(0.05), col = "gray", lty = 2)
abline(h = -log10(fdrcut), col = 2, lty = 2)
points(resNS$beta.1 * 100, -log10(resNS$pvalue), pch = 21, col = ifelse(resNS$diffage, "black", cell2col(NS$celltype)), bg = ifelse(resNS$diffage, cell2col(NS$celltype), "white"))
box()

i <- which.min(resNS$pvalue)
text(resNS$beta.1[i]*100, -log10(resNS$pvalue[i]), rownames(resNS)[i], pos = 2)

axis(1)
axis(2)

title("Schizophrenia", adj = 0, line = 0.5)
title(xlab = bquote(beta ~ "coefficient"))
title(ylab = "-Log10 p-value", line = 1.6)

mtext("b", 3, adj = -0.2, line = 1.7, font = 2)

# schizophrenia T-cells
plot.new()
plot.window(xlim = range(resTS$beta.1 * 100), ylim = c(0,6))
abline(h = -log10(0.05), col = "gray", lty = 2)
abline(h = -log10(fdrcut), col = 2, lty = 2)
points(resTS$beta.1 * 100, -log10(resTS$pvalue), pch = 21, col = ifelse(resTS$diffage, "black", cell2col("tcd4")), bg = ifelse(resTS$diffage, cell2col("tcd4"), "white"))
box()

axis(1)
axis(2)

title(xlab = bquote(beta ~ "coefficient"))
title(ylab = "-Log10 p-value", line = 1.6)

legend("bottomright", c("Cell-DEA:", "yes", "no"), text.width = strwidth(c("Cell-DEA:", "yes", "no")), x.intersp = 0.5, pch = c(NA,19,21), col = c(NA, "gray40", "gray20"), bty = "n", inset = c(0, 0.99), xpd = TRUE, horiz = TRUE)


# BMI neutrophils
fdrcut <- mean(c(max(resNB$pvalue[p.adjust(resNB$pvalue, "fdr") < 0.05]), min(resNB$pvalue[p.adjust(resNB$pvalue, "fdr") > 0.05])))

plot.new()
plot.window(xlim = c(-0.4,0.4), ylim = c(0,5.5))
abline(h = -log10(0.05), col = "gray", lty = 2)
abline(h = -log10(fdrcut), col = 2, lty = 2)
points(resNB$beta.1 * 100, -log10(resNB$pvalue), pch = 21, col = ifelse(resNB$diffage, "black", cell2col(NB$celltype)), bg = ifelse(resNB$diffage, cell2col(NB$celltype), "white"))
box()

i <- order(resNB$pvalue)[1]
text(resNB$beta.1[i]*100, -log10(resNB$pvalue[i]), rownames(resNB)[i], adj = c(0.8,-0.8))

i <- order(resNB$pvalue)[2]
text(resNB$beta.1[i]*100, -log10(resNB$pvalue[i]), rownames(resNB)[i], pos = 2)

axis(1)
axis(2)

title("BMI", adj = 0, line = 0.5)
title(xlab = bquote(beta ~ "coefficient"))
title(ylab = "-Log10 p-value", line = 1.6)

mtext("c", 3, adj = -0.2, line = 1.7, font = 2)


# BMI T-cells
fdrcut <- mean(c(max(resTB$pvalue[p.adjust(resTB$pvalue, "fdr") < 0.05]), min(resTB$pvalue[p.adjust(resTB$pvalue, "fdr") > 0.05])))

plot.new()
plot.window(xlim = c(-0.4,0.4), ylim = c(0,5.5))
abline(h = -log10(0.05), col = "gray", lty = 2)
abline(h = -log10(fdrcut), col = 2, lty = 2)
points(resTB$beta.1 * 100, -log10(resTB$pvalue), pch = 21, col = ifelse(resTB$diffage, "black", cell2col("tcd4")), bg = ifelse(resTB$diffage, cell2col("tcd4"), "white"))
box()

i <- order(resTB$pvalue)[1]
text(resTB$beta.1[i]*100, -log10(resTB$pvalue[i]), rownames(resTB)[i], adj = c(0,-0.8))

i <- order(resTB$pvalue)[2]
text(resTB$beta.1[i]*100, -log10(resTB$pvalue[i]), rownames(resTB)[i], adj = c(0,-0.6))

axis(1)
axis(2)

title(xlab = bquote(beta ~ "coefficient"))
title(ylab = "-Log10 p-value", line = 1.6)

legend("bottomright", c("Cell-DEA:", "yes", "no"), text.width = strwidth(c("Cell-DEA:", "yes", "no")), x.intersp = 0.5, pch = c(NA,19,21), col = c(NA, "gray40", "gray20"), bty = "n", inset = c(0, 0.99), xpd = TRUE, horiz = TRUE)


# smoking neutrophils
fdrcut <- mean(c(max(resNT$pvalue[p.adjust(resNT$pvalue, "fdr") < 0.05]), min(resNT$pvalue[p.adjust(resNT$pvalue, "fdr") > 0.05])))

plot.new()
plot.window(xlim = c(-25,25), ylim = c(0,18))
abline(h = -log10(0.05), col = "gray", lty = 2)
abline(h = -log10(fdrcut), col = 2, lty = 2)
points(resNT$beta.1 * 100, -log10(resNT$pvalue), pch = 21, col = ifelse(resNT$diffage, "black", cell2col(NT$celltype)), bg = ifelse(resNT$diffage, cell2col(NT$celltype), "white"))
box()

i <- order(resNT$pvalue)[1]
text(resNT$beta.1[i]*100, -log10(resNT$pvalue[i]), rownames(resNT)[i], adj = c(-0.1,0.3))

i <- order(resNT$pvalue)[2]
text(resNT$beta.1[i]*100, -log10(resNT$pvalue[i]), rownames(resNT)[i], adj = c(0,-0.8))

i <- order(resNT$pvalue)[3]
text(resNT$beta.1[i]*100, -log10(resNT$pvalue[i]), rownames(resNT)[i], adj = c(0,1.3))

axis(1)
axis(2)

title("Smoking", adj = 0, line = 0.5)
title(xlab = bquote(beta ~ "coefficient"))
title(ylab = "-Log10 p-value", line = 1.6)

mtext("d", 3, adj = -0.2, line = 1.7, font = 2)


# smoking T-cells
fdrcut <- mean(c(max(resTT$pvalue[p.adjust(resTT$pvalue, "fdr") < 0.05]), min(resTT$pvalue[p.adjust(resTT$pvalue, "fdr") > 0.05])))

plot.new()
plot.window(xlim = c(-25,25), ylim = c(0,18))
abline(h = -log10(0.05), col = "gray", lty = 2)
abline(h = -log10(fdrcut), col = 2, lty = 2)
points(resTT$beta.1 * 100, -log10(resTT$pvalue), pch = 21, col = ifelse(resTT$diffage, "black", cell2col("tcd4")), bg = ifelse(resTT$diffage, cell2col("tcd4"), "white"))
box()

i <- order(resTT$pvalue)[1]
text(resTT$beta.1[i]*100, -log10(resTT$pvalue[i]), rownames(resTT)[i], adj = c(-0.05,-0.7))

i <- order(resTT$pvalue)[2]
text(resTT$beta.1[i]*100, -log10(resTT$pvalue[i]), rownames(resTT)[i], adj = c(-0.05,-0.7))

i <- order(resTT$pvalue)[3]
text(resTT$beta.1[i]*100, -log10(resTT$pvalue[i]), rownames(resTT)[i], adj = c(-0.05,-0.7))

axis(1)
axis(2)

title(xlab = bquote(beta ~ "coefficient"))
title(ylab = "-Log10 p-value", line = 1.6)

legend("bottomright", c("Cell-DEA:", "yes", "no"), text.width = strwidth(c("Cell-DEA:", "yes", "no")), x.intersp = 0.5, pch = c(NA,19,21), col = c(NA, "gray40", "gray20"), bty = "n", inset = c(0, 0.99), xpd = TRUE, horiz = TRUE)


# neutrophil overlaps
tabNS <- table(resNS$pvalue <= 0.05 & sign(resNS$beta.1) == sign(ewasS$diff), resNS$diffage)
tabNB <- table(resNB$pvalue <= 0.05 & sign(resNB$beta.1) == sign(ewasB$diff), resNB$diffage)
tabNT <- table(resNT$pvalue <= 0.05 & sign(resNT$beta.1) == sign(ewasT$diff), resNT$diffage)

repN <- c((tabNS[2,]/colSums(tabNS))[2:1], (tabNB[2,]/colSums(tabNB))[2:1], (tabNT[2,]/colSums(tabNT))[2:1]) * 100

plot.new()
plot.window(xlim = c(0.5,8.6), ylim = c(0,60), yaxs = "i")
abline(h = seq(0,60,10), col = "gray80", lty = 3)
rect(c(1,2,4,5,7,8)-0.25, 0, c(1,2,4,5,7,8)+0.25, repN, col = rep(c(cell2col("neu"), "white"), 3), border = c("black", cell2col("neu")))

title(ylab = "Replicated hits (%)", line = 1.6)
axis(1, c(1.5,4.5,7.5), c("Schizophrenia", "BMI", "Smoking"))
axis(2)
box(bty = "l")

mtext("e", 3, adj = -0.2, line = 1.7, font = 2)


# tcell overlaps
tabTS <- table(resTS$pvalue <= 0.05 & sign(resTS$beta.1) == sign(ewasS$diff), resTS$diffage)
tabTB <- table(resTB$pvalue <= 0.05 & sign(resTB$beta.1) == sign(ewasB$diff), resTB$diffage)
tabTT <- table(resTT$pvalue <= 0.05 & sign(resTT$beta.1) == sign(ewasT$diff), resTT$diffage)

repT <- c((tabTS[2,]/colSums(tabTS))[2:1], (tabTB[2,]/colSums(tabTB))[2:1], (tabTT[2,]/colSums(tabTT))[2:1]) * 100

plot.new()
plot.window(xlim = c(0.5,8.6), ylim = c(0,60), yaxs = "i")
abline(h = seq(0,60,10), col = "gray80", lty = 3)
rect(c(1,2,4,5,7,8)-0.25, 0, c(1,2,4,5,7,8)+0.25, repT, col = rep(c(cell2col("tcd4"), "white"), 3), border = c("black", cell2col("tcd4")))

title(ylab = "Replicated hits (%)", line = 1.6)
axis(1, c(1.5,4.5,7.5), c("Schizophrenia", "BMI", "Smoking"))
axis(2)
box(bty = "l")

legend("bottomright", c("Cell-DEA:", "yes", "no"), text.width = strwidth(c("Cell-DEA:", "yes", "no")), x.intersp = 0.5, fill = c(NA, "gray40", "white"), border = c(NA, "gray40", "gray20"), bty = "n", inset = c(0, 0.99), xpd = TRUE, horiz = TRUE)

invisible(dev.off())
