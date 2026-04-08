#--- init ----------------------------------------------------------------------

library(annmatrix)
library(kklibrary)
suppressPackageStartupMessages(library(inops))

dfage <- readRDS("data/diffage_cgs_hannon.rds")
ewas  <- read.csv("data/EWAS/ewas_alz.csv", skip = 2)


#--- clean ---------------------------------------------------------------------

rownames(ewas) <- ewas[,1]
ewas <- ewas[,c(6,7)]
colnames(ewas) <- c("diff", "pvalue")

cgs   <- intersect(rownames(ewas), rownames(dfage))
ewas  <- ewas[cgs,]
dfage <- dfage[cgs,]

# select only diff age
cgs   <- dfage$category == "cell-dea"
ewas  <- ewas[cgs,]
dfage <- dfage[cgs,]


#--- plot ----------------------------------------------------------------------

# NOTE: manual colors
cell2col <- function(x) {
  cols <- as.character(x)
  cols[] <- "black"
  cols[x == "neun+"]   <- "gray60"
  cols[x == "sox10+"]  <- "#8D4B08"
  cols[x == "double-"] <- "black"
  cols
}

myplot <- function(age, diff, col, ...) {
  plot(age, diff, col = col, pch = 19, cex = 0.6, xlab = "", ylab = "", axes = FALSE, ...)
  abline(h = 0, col = "gray", lty = 2)
  abline(v = 0, col = "gray", lty = 2)
  title(xlab = bquote(Delta ~ "Modification (%/year)"))
  title(ylab = bquote("Plaque burden" ~ beta ~ "coefficient"), line = 1.25)
  box()

  fish1 <- fisher.test(table(age < 0, diff > 0))
  fish2 <- fisher.test(table(age > 0, diff > 0))
  fish  <- if(fish1$estimate > fish2$estimate) fish1 else fish2
  fish$estimate <- ifelse(fish$estimate > 100, round(fish$estimate), ifelse(fish$estimate > 10, round(fish$estimate, 1), round(fish$estimate, 2)))
  fish$p.value  <- ifelse(fish$p.value < 0.001, formatC(fish$p.value, format = "e", digits = 2), signif(fish$p.value, 2))
  legend("topleft", paste(c("OR =", "p ="), c(fish$estimate, fish$p.value)), text.col = "red2", bty = "n", inset = c(-0.1,-0.03))
}

pdf("fig_8_c.pdf", width = 15/2.54, height = 4/2.54, pointsize = 8)

par(mfrow = c(1,4), mar = c(3,3,3,0), oma = c(0,0,0,3), tck = -0.02, mgp = c(2,0.5,0), las = 1)

myplot(dfage$agingrate_neunpos, ewas$diff, cell2col("neun+"), xlim = c(-0.25,0.25), ylim = c(-8,11))
axis(1, at = c(-0.2, 0, 0.2))
axis(2, at = c(-6,0,6), las = 2)
mtext("c", 3, adj = -0.18, line = 1.7, font = 2)

myplot(dfage$agingrate_sox10pos, ewas$diff, cell2col("sox10+"), xlim = c(-0.25,0.25), ylim = c(-8,11))
axis(1, at = c(-0.2, 0, 0.2))
axis(2, at = c(-6,0,6), las = 2)

myplot(dfage$agingrate_doubleneg, ewas$diff, cell2col("double-"), xlim = c(-1,1), ylim = c(-8,11))
axis(1, at = c(-0.8, 0, 0.8))
axis(2, at = c(-6,0,6), las = 2)

legend("bottom", legend = c("NeuN+", "SOX10+", "Double-"), text.width = strwidth(c("NeuN+", "SOX10+", "Double-")), col = cell2col(c("neun+", "sox10+", "double-")),
       pch = 19, bty = 'n', inset = c(0, 1.025), xpd = NA, horiz = TRUE)

invisible(dev.off())
