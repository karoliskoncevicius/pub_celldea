#--- init ----------------------------------------------------------------------

library(stringr)
library(annmatrix)
library(basetheme)
library(kklibrary)

X <- readRDS("data/blueprint/3_annmatrix.rds")


#--- plot ----------------------------------------------------------------------

# so that celltypes do not overlap each other on the plot
X$celltype <- droplevels(factor(X$celltype, levels = c("neu", "mono", "bcell", "bn", "bm", "nk", "tcell", "tcd4n", "tcd4m", "tcd8n", "tcd8m")))
X <- X[, sample(1:ncol(X))]

# NOTE: manually selected cytosines
cgs <- c("cg18263572", "cg10501210", "cg00751770", "cg06893560", "cg15473346")


pdf("fig_2_a.pdf", width = 15/2.54, height = 4/2.54, pointsize = 8)

par(mfrow = c(1,5), mar = c(4,0,4,1), oma = c(0,4,0,0), tck = -0.02, mgp = c(2,0.5,0))

for(cg in cgs) {
  if(cg  %in% rownames(X)) {
    plot(X$age, X[cg,] * 100, col = cell2col(X$celltype), pch = 19, cex = 0.8, xlim = c(19,81), ylim = c(0,100), xlab = "", yaxt = "n", yaxs = "i")
  } else {
    plot.new()
  }
  title(cg, line = 0.5)
  if(cg == cgs[1]) {
    mtext("a", 3, line = 2, font = 2, adj = -0.3)
    title(ylab = "Modification (%)", line = 2.5, outer = TRUE)
    axis(2, las = 2)
  }
  if(cg == cgs[3]) {
    title(xlab = "Age (years)")
  }
}

legend("bottom", legend = c("Neu", "Mono", "TCD4n"), text.width = strwidth(c("Neu", "Mono", "TCD4n")),
       col = cell2col(c("neu", "mono", "tcd4n")), pch = 19, bty = 'n', inset = c(0, 1.15), xpd = NA, horiz = TRUE)

invisible(dev.off())
