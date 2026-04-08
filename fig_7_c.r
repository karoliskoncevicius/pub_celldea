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

# function
myplot <- function(age, diff, col, ...) {
  plot(age, diff, col = adjustcolor(col, 0.25), pch = 19, cex = 0.6, xlab = "", ylab = "", axes = FALSE, ...)
  abline(h = 0, col = "gray", lty = 2)
  abline(v = 0, col = "gray", lty = 2)
  title(xlab = bquote(Delta ~ "Modification (%/year)"))
  title(ylab = bquote("SCZ" ~ beta ~ "coefficient"), line = 1.25)
  axis(1, at = c(-0.6, 0, 0.6))
  axis(2, at = c(-2,0,2), las = 2)
  box()

  fish1 <- fisher.test(table(age < 0, diff > 0))
  fish2 <- fisher.test(table(age > 0, diff > 0))
  fish  <- if(fish1$estimate > fish2$estimate) fish1 else fish2
  fish$estimate <- ifelse(fish$estimate > 100, round(fish$estimate), ifelse(fish$estimate > 10, round(fish$estimate, 1), round(fish$estimate, 2)))
  fish$p.value  <- ifelse(fish$p.value == 0, "< 1e-100", ifelse(fish$p.value < 0.001, formatC(fish$p.value, format = "e", digits = 1), signif(fish$p.value, 2)))
  legend("topleft", paste(c("OR =", "p ="), c(fish$estimate, fish$p.value)), text.col = "red2", bty = "n", inset = c(-0.1,-0.03))
}

pdf("fig_7_c.pdf", width = 14/2.54, height = 4/2.54, pointsize = 8)

par(mfrow = c(1,4), mar = c(3,3,3,0), oma = c(0,0,0,1), tck = -0.02, mgp = c(2,0.5,0), las = 1)

inds <- df$category == "cell-dea" & p.adjust(fits[[1]]$pvalue, "fdr") <= 0.05
myplot(df$ra_agingrate_tcd4n[inds], fits[[1]]$beta.1[inds]*100, cell2col("tcd4n"), xlim = c(-0.9,0.9), ylim = c(-3,3))
title("EpiDISH-7", line = 0.5)
mtext("c", 3, adj = -0.2, line = 1.7, font = 2)

inds <- df$category == "cell-dea" & p.adjust(fits[[3]]$pvalue, "fdr") <= 0.05
myplot(df$ra_agingrate_tcd4n[inds], fits[[3]]$beta.1[inds]*100, cell2col("tcd4n"), xlim = c(-0.9,0.9), ylim = c(-3,3))
title("EpiDISH-12", line = 0.5)

inds <- df$category == "cell-dea" & p.adjust(fits[[5]]$pvalue, "fdr") <= 0.05
myplot(df$ra_agingrate_tcd4n[inds], fits[[5]]$beta.1[inds]*100, cell2col("tcd4n"), xlim = c(-0.9,0.9), ylim = c(-3,3))
title("DEA-free", line = 0.5, xpd = NA)

inds <- df$category == "cell-dea" & p.adjust(fits[[6]]$pvalue, "fdr") <= 0.05
myplot(df$ra_agingrate_tcd4n[inds], fits[[6]]$beta.1[inds]*100, cell2col("tcd4n"), xlim = c(-0.9,0.9), ylim = c(-3,3))
title("DEA-free (celltype x age)", line = 0.5, xpd = NA)

invisible(dev.off())
