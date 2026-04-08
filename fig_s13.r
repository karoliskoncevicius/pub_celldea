#--- init ----------------------------------------------------------------------

library(annmatrix)
library(kklibrary)

fits <- readRDS("data/scz_ewas_correction.rds")
df   <- readRDS("data/diffage_cgs.rds")


#--- merge ---------------------------------------------------------------------

cgs <- intersect(rownames(fits[[1]]), rownames(df))

df <- df[cgs,]

for(i in 1:length(fits)) {
  fits[[i]] <- fits[[i]][cgs,]
}

stopifnot(all.equal(rownames(fits[[1]]), rownames(df)))


#--- plot ----------------------------------------------------------------------

pdf("fig_s13.pdf", width = 4/2.54, height = 4/2.54, pointsize = 8)

par(mar = c(0,1,0,0), tck = -0.02, mgp = c(2,0.5,0), las = 1)

i1 <- p.adjust(fits[[3]]$pvalue, "fdr") <= 0.05
i2 <- p.adjust(fits[[6]]$pvalue, "fdr") <= 0.05

v1 <- i1 & !i2
v2 <- !i1 & i2
v3 <- i1 & i2


symbols(c(0.65,1.35), c(1,1), c(0.7,0.7), inches = FALSE, xlim = c(0,2.1), lwd = 2, ann = FALSE, axes = FALSE)
text(c(0.5,1.5), c(1.35,1.355), c("EpiDISH", "DEA-free"), font = 2)
text(c(0.4,1,1.6), c(1.05,1.05,1.05), c(sum(v1),sum(v3),sum(v2)))
text(c(0.35,1,1.65), c(0.925,0.925,0.925), paste0(round(c(mean(df$category[v1] == "cell-dea"), mean(df$category[v3] == "cell-dea"), mean(df$category[v2] == "cell-dea"))*100, 1), " %"))

invisible(dev.off())
