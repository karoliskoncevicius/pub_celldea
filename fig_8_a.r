#--- init ----------------------------------------------------------------------

library(annmatrix)
library(basetheme)
library(kklibrary)

X <- readRDS("data/hannon/3_annmatrix.rds")


#--- plot ----------------------------------------------------------------------

# so that celltypes do not overlap each other on the plot
X$celltype[X$celltype == "neunpos"]   <- "NeuN+"
X$celltype[X$celltype == "sox10pos"]  <- "SOX10+"
X$celltype[X$celltype == "doubleneg"] <- "Double-"
X$celltype <- droplevels(factor(X$celltype, levels = c("NeuN+", "SOX10+", "Double-")))
X <- X[, sample(1:ncol(X))]

# NOTE: manually selected cytosines
cgs <- c("cg12935715", "cg25209674", "cg24724428", "cg27344602", "cg18241597")

# NOTE: manual colors
cell2col <- function(x) {
  cols <- as.character(x)
  cols[] <- "black"
  cols[x == "NeuN+"]   <- "gray60"
  cols[x == "SOX10+"]  <- "#8D4B08"
  cols[x == "Double-"] <- "black"
  cols
}


pdf("fig_8_a.pdf", width = 15/2.54, height = 4/2.54, pointsize = 8)

par(mfrow = c(1,5), mar = c(4,0,4,1), oma = c(0,4,0,0), tck = -0.02, mgp = c(2,0.5,0))

for(cg in cgs) {
  if(cg  %in% rownames(X)) {
    plot(X$age, X[cg,] * 100, col = cell2col(X$celltype), pch = 19, cex = 0.8, xlim = c(17,100), ylim = c(0,100), xlab = "", yaxt = "n", yaxs = "i")
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

legend("bottom", legend = levels(X$celltype), text.width = strwidth(levels(X$celltype)), col = cell2col(levels(X$celltype)), pch = 19, bty = 'n', inset = c(0, 1.15), xpd = NA, horiz = TRUE)

invisible(dev.off())
