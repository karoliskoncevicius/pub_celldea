#--- init ----------------------------------------------------------------------

library(stringr)
library(annmatrix)
library(kklibrary)
library(matrixTests)
suppressPackageStartupMessages(library(inops))

res <- readRDS("data/ewasoverlaps_blueprint.rds")


#--- plot ----------------------------------------------------------------------

res$pmid  <- factor(res$pmid, levels = unique(res$pmid))
res$trait <- str_to_sentence(res$trait)

format_pvals <- function(p) {
  if(p < 1*10^-100) {
    "< 1e-100"
  } else if(p < 0.0001) {
    sprintf("%0.1e", p)
  } else {
    signif(p, 2)
  }
}

res$odds_rest[res$odds_rest < 1/2^8] <- 1/2^8


pdf("fig_s3.pdf", width = 12/2.54, height = nlevels(res$pmid)/4.5/2.54, pointsize = 8)

layout(matrix(1:3, nrow = 1), widths=c(0.35,0.4,0.25))

par(mar=c(3,0,2,0), tck=-0.02, mgp=c(2,0.5,0))
plot.new()
plot.window(xlim=c(-2.5, 1), ylim=c(nlevels(res$pmid)+0.5, 0.5), xaxs="i", yaxs="i")

strps <- seq(1, nlevels(res$pmid)+1, 2)
rect(-2.35, strps-0.5, 100, strps+0.5, col="gray92", border=NA)

for(i in 1:nlevels(res$pmid)) {
  ind <- which(res$pmid == levels(res$pmid)[i])[1]
  text(1, i+0.15, bquote(bold(.(res$trait[ind])) ~ "("*.(res$author[ind]) ~ .(res$year[ind])*")"), pos = 2)
}

par(mar=c(3,0,2,0), tck=-0.02, mgp=c(2,0.5,0))
plot.new()
plot.window(xlim = c(-9,6.2), ylim=c(nlevels(res$pmid)+0.5, 0.5), xaxs="i", yaxs="i")

strps <- seq(1, nlevels(res$pmid)+1, 2)
rect(-100, strps-0.5, 100, strps+0.5, col = "gray92", border = NA)
abline(v = 0, lty = 3, col = "black")

points(log2(res$odds_rest), res$pmid, pch = ifelse(res$pval_rest <= 0.05, 19, 1), col = "darkgray")
points(log2(res$odds_cell), res$pmid, pch = ifelse(res$pval_cell <= 0.05, 19, 1), col = "blue2")
points(log2(res$odds_diff), res$pmid, pch = ifelse(res$pval_diff <= 0.05, 19, 1), col = "red2")

axis(1, at = seq(-6,6,2))
axis(1, at = -8, labels = c("< -6"))
mtext("Log2 odds ratio", 1, at = 0, line = 1.5, cex = 0.6)


par(mar = c(3,0,2,1), tck = -0.02, mgp = c(2,0.5,0))
plot.new()
plot.window(xlim = c(0.5, 3.5), ylim = c(nlevels(res$pmid)+0.5, 0.5), xaxs = "i", yaxs = "i")

strps <- seq(1, nlevels(res$pmid)+1, 2)
rect(-100, strps-0.5, 100, strps+0.5, col = "gray92", border = NA)

text(1:3, 0.25, c("Percent", "Odds ratio", "p-value"), pos = 3, xpd = TRUE, font = 2)
text(1, 1:nrow(res), paste0(sprintf("%0.1f", res$perc_diff[!duplicated(res$pmid)]), " %"))
text(2, 1:nrow(res), sprintf("%0.1f", res$odds_diff[!duplicated(res$pmid)]))
text(3, 1:nrow(res), mapply(format_pvals, res$pval_diff[!duplicated(res$pmid)]))

invisible(dev.off())
