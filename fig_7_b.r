#--- init ----------------------------------------------------------------------

library(readxl)
library(annmatrix)
library(kklibrary)

fits <- readRDS("data/scz_ewas_correction.rds")
df   <- readRDS("data/diffage_cgs.rds")

ewas <- as.data.frame(read_excel("data/EWAS/elife-58430-supp1-v2.xlsx", sheet = 3, skip = 3))


#--- clean ---------------------------------------------------------------------

rownames(ewas) <- ewas[,1]
ewas <- ewas[,c(7,9)]
colnames(ewas) <- c("diff", "pvalue")


#--- merge ---------------------------------------------------------------------

cgs <- intersect(rownames(fits[[1]]), rownames(df))

df <- df[cgs,]

for(i in 1:length(fits)) {
  fits[[i]] <- fits[[i]][cgs,]
}

stopifnot(all.equal(rownames(fits[[1]]), rownames(df)))


#--- overlaps ------------------------------------------------------------------

sigs <- list()
for(f in names(fits)) {
  sigs[[paste0(f, "_fdr")]] <- p.adjust(fits[[f]]$pvalue, "fdr") <= 0.05
}

res <- data.frame(id = names(sigs))
for(s in names(sigs)) {
  is_sig <- sigs[[s]]
  res$n[res$id == s] <- sum(is_sig)
  # static
  fish <- fisher.test(table(is_sig, df$category == "age-static"))
  res$odds_cell[res$id == s] <- fish$estimate
  res$pval_cell[res$id == s] <- fish$p.value
  res$perc_cell[res$id == s] <- mean((df$category == "age-static")[is_sig]) * 100
  # diffage
  fish <- fisher.test(table(is_sig, df$category == "cell-dea"))
  res$odds_diff[res$id == s] <- fish$estimate
  res$pval_diff[res$id == s] <- fish$p.value
  res$perc_diff[res$id == s] <- mean((df$category == "cell-dea")[is_sig]) * 100
  # rest
  fish <- fisher.test(table(is_sig, df$category == "rest"))
  res$odds_rest[res$id == s] <- fish$estimate
  res$pval_rest[res$id == s] <- fish$p.value
  res$perc_rest[res$id == s] <- mean((df$category == "rest")[is_sig]) * 100
}

res$reference <- NA
res$reference[startsWith(res$id, "epi_6ct")]  <- "EpiDISH-7"
res$reference[startsWith(res$id, "epi_12ct")] <- "EpiDISH-12"
res$reference[startsWith(res$id, "new_12ct")] <- "DEA-free"


#--- plot ----------------------------------------------------------------------

format_pvals <- function(p) {
  if(p < 1*10^-100) {
    "< 1e-100"
  } else if(p < 0.0001) {
    sprintf("%0.1e", p)
  } else {
    signif(p, 2)
  }
}

pdf("fig_7_b.pdf", width = 14/2.54, height = nrow(res)/2/2.54, pointsize = 8)

layout(matrix(1:3, nrow = 1), widths=c(0.30,0.45,0.25))

par(mar = c(3,0,3,0), tck = -0.02, mgp = c(2,0.5,0))
plot.new()
plot.window(xlim=c(-2.5, 1.5), ylim=c(nrow(res)+0.5, 0.5), xaxs="i", yaxs="i")

strps <- seq(1, nrow(res)+1, 2)
rect(-2.35, strps-0.5, 100, strps+0.5, col="gray92", border=NA)

text(c(-0.7,0.4,1.5), -0.25, c("Reference", "Model", "# CpGs"), pos = 2, xpd = TRUE, font = 2)
text(-0.7, (1:nrow(res))+0.15, res$reference, pos = 2)
text(0.4, (1:nrow(res))+0.15, ifelse(grepl("age", res$id), "CT x age", "CT + age"), pos = 2)
text(1.5, (1:nrow(res))+0.15, res$n, pos = 2)

mtext("b", 3, adj = 0.05, line = 1, font = 2)

par(mar=c(3,0,3,0), tck=-0.02, mgp=c(2,0.5,0))
plot.new()
plot.window(xlim = c(-4.5,4.5), ylim=c(nrow(res)+0.5, 0.5), xaxs="i", yaxs="i")

strps <- seq(1, nrow(res)+1, 2)
rect(-100, strps-0.5, 100, strps+0.5, col = "gray92", border = NA)
abline(v = 0, lty = 3, col = "black")

points(log2(res$odds_diff), 1:nrow(res), pch = ifelse(res$pval_diff <= 0.05, 19, 1), col = "red2")
points(log2(res$odds_cell), 1:nrow(res), pch = ifelse(res$pval_cell <= 0.05, 19, 1), col = "blue2")
points(log2(res$odds_rest), 1:nrow(res), pch = ifelse(res$pval_rest <= 0.05, 19, 1), col = "darkgray")

axis(1, at = seq(-4,4,2))
mtext("Log2 odds ratio", 1, at = 0, line = 1.5, cex = 0.6)


par(mar = c(3,0,3,1), tck = -0.02, mgp = c(2,0.5,0))
plot.new()
plot.window(xlim = c(0.5, 3.5), ylim = c(nrow(res)+0.5, 0.5), xaxs = "i", yaxs = "i")

strps <- seq(1, nrow(res)+1, 2)
rect(-100, strps-0.5, 100, strps+0.5, col = "gray92", border = NA)

text(1:3, 0.25, c("Percent", "Odds ratio", "p-value"), pos = 3, xpd = TRUE, font = 2)
text(1, 1:nrow(res), paste0(sprintf("%0.1f", res$perc_diff), " %"))
text(2, 1:nrow(res), sprintf("%0.2f", res$odds_diff))
text(3, 1:nrow(res), mapply(format_pvals, res$pval_diff))

invisible(dev.off())
