#--- init ----------------------------------------------------------------------

library(stringr)
library(annmatrix)
library(kklibrary)
library(matrixStats)
suppressPackageStartupMessages(library(inops))

dfage <- readRDS("data/diffage_cgs.rds")
ewas  <- readRDS("data/ewascatalog_blood.rds")
res   <- readRDS("data/ewasoverlaps.rds")


#--- filter --------------------------------------------------------------------

# select top 20
res <- res[res$pmid %in% unique(res$pmid)[1:20],]

# match ewas
ewas <- ewas[ewas$studyid %in% res$study,]

# turn off because of minue negative infinity overlaps for all cells that are not significant
ewas$beta[ewas$studyid == "31073081_Imboden-M_fev1_meta-analysis"] <- NA


#--- overlaps ------------------------------------------------------------------

is_diff <- dfage$category == "cell-dea"

myfish <- function(x, y) {
  if(length(x) > 1) {
    fish <- fisher.test(factor(sign(x), levels = c(-1,1)), factor(sign(y), levels = c(-1,1)))
    fish$percent <- mean(sign(x) == sign(y)) * 100
  } else {
    fish <- list(estimate = NA, p.value = NA, percent = NA)
  }
  fish
}

# add new columns to res
for(s in res$study) {

  e <- ewas[ewas$studyid == s,]
  e <- e[!duplicated(e$cpg),]
  rownames(e) <- e$cpg

  if(!all(is.na(e$beta))) {

    cgs  <- intersect(rownames(e), rownames(dfage)[is_diff])
    inds <- res$study == s

    fish <- myfish(e[cgs,]$beta, dfage[cgs,]$ra_agingrate_tcd4n)
    res$diff_odds_tcd4n[inds] <- fish$estimate
    res$diff_pval_tcd4n[inds] <- fish$p.value
    res$diff_perc_tcd4n[inds] <- ifelse(fish$percent < 50, 100 - fish$percent, fish$percent)

    fish <- myfish(e[cgs,]$beta, dfage[cgs,]$ra_agingrate_tcd4m)
    res$diff_odds_tcd4m[inds] <- fish$estimate
    res$diff_pval_tcd4m[inds] <- fish$p.value
    res$diff_perc_tcd4m[inds] <- ifelse(fish$percent < 50, 100 - fish$percent, fish$percent)

    fish <- myfish(e[cgs,]$beta, dfage[cgs,"ra_agingrate_mono"])
    res$diff_odds_mono[inds] <- fish$estimate
    res$diff_pval_mono[inds] <- fish$p.value
    res$diff_perc_mono[inds] <- ifelse(fish$percent < 50, 100 - fish$percent, fish$percent)

    fish <- myfish(e[cgs,]$beta, dfage[cgs,]$bp_agingrate_neu)
    res$diff_odds_neu[inds] <- fish$estimate
    res$diff_pval_neu[inds] <- fish$p.value
    res$diff_perc_neu[inds] <- ifelse(fish$percent < 50, 100 - fish$percent, fish$percent)

    i  <- which.max(abs(log2(c(res$diff_odds_tcd4n[inds], res$diff_odds_tcd4m[inds], res$diff_odds_mono[inds], res$diff_odds_neu[inds]))))
    cs <- rank(abs(log2(c(res$diff_odds_tcd4n[inds], res$diff_odds_tcd4m[inds], res$diff_odds_mono[inds], res$diff_odds_neu[inds]))))
    res$diff_max_cell[inds] <- paste(c("tcd4n", "tcd4m", "mono", "neu")[cs == max(cs)], collapse = "/")
    res$diff_max_odds[inds] <- c(res$diff_odds_tcd4n[inds], res$diff_odds_tcd4m[inds], res$diff_odds_mono[inds], res$diff_odds_neu[inds])[i]
    res$diff_max_pval[inds] <- c(res$diff_pval_tcd4n[inds], res$diff_pval_tcd4m[inds], res$diff_pval_mono[inds], res$diff_pval_neu[inds])[i]
    res$diff_max_perc[inds] <- c(res$diff_perc_tcd4n[inds], res$diff_perc_tcd4m[inds], res$diff_perc_mono[inds], res$diff_perc_neu[inds])[i]

    cgs <- intersect(rownames(e), rownames(dfage)[!is_diff])

    fish <- myfish(e[cgs,]$beta, dfage[cgs,]$ra_agingrate_tcd4n)
    res$stat_odds_tcd4n[res$study == s] <- fish$estimate
    res$stat_pval_tcd4n[res$study == s] <- fish$p.value

    fish <- myfish(e[cgs,]$beta, dfage[cgs,]$ra_agingrate_tcd4m)
    res$stat_odds_tcd4m[res$study == s] <- fish$estimate
    res$stat_pval_tcd4m[res$study == s] <- fish$p.value

    fish <- myfish(e[cgs,]$beta, dfage[cgs,]$ra_agingrate_mono)
    res$stat_odds_mono[res$study == s] <- fish$estimate
    res$stat_pval_mono[res$study == s] <- fish$p.value

    fish <- myfish(e[cgs,]$beta, dfage[cgs,]$bp_agingrate_neu)
    res$stat_odds_neu[res$study == s] <- fish$estimate
    res$stat_pval_neu[res$study == s] <- fish$p.value
  }
}

# select only those with non-NA overlaps (those that provided beta values)
res <- res[!is.na(res$diff_odds_tcd4n),]


#--- plot ----------------------------------------------------------------------

res$pmid <- factor(res$pmid, levels = unique(res$pmid))
res$trait <- str_to_sentence(res$trait)

format_pvals <- function(p) {
  if(p < 1*10^-100) {
    "< 1e-100"
  } else if(p < 0.001) {
    sprintf("%0.1e", p)
  } else {
    signif(p, 2)
  }
}

format_perc <- function(p) {
  if(is.na(p)) {
    ""
  } else {
    paste0(round(p), " %")
  }
}

format_odds <- function(x) {
  x <- 2^abs(-log2(x))
  if(x > 100) {
    round(x)
  } else {
    round(x, 1)
  }
}

cap_or <- function(x) {
  x[which(log2(x) > 11)]  <- 2^13
  x[which(log2(x) < -11)] <- 2^-13
  x
}


pdf("fig_4_c.pdf", width = 14/2.54, height = nlevels(res$pmid)/3.25/2.54, pointsize = 8)

layout(matrix(1:3, nrow = 1), widths=c(0.30,0.45,0.25))

par(mar=c(3,0,2,0), tck=-0.02, mgp=c(2,0.5,0))
plot.new()
plot.window(xlim=c(-2.5, 1), ylim=c(nlevels(res$pmid)+0.5, 0.5), xaxs="i", yaxs="i")

strps <- seq(1, nlevels(res$pmid)+1, 2)
rect(-2.35, strps-0.5, 100, strps+0.5, col="gray92", border=NA)

for(i in 1:nlevels(res$pmid)) {
  ind <- which(res$pmid == levels(res$pmid)[i])[1]
  text(1, i+0.15, bquote(bold(.(res$trait[ind])) ~ "("*.(res$author[ind]) ~ .(res$year[ind])*")"), pos = 2)
}

mtext("c", 3, adj = 0.03, line = 1, font = 2)

par(mar=c(3,0,2,0), tck=-0.02, mgp=c(2,0.5,0))
plot.new()
plot.window(xlim = c(-11,15.5), ylim=c(nlevels(res$pmid)+0.5, 0.5), xaxs="i", yaxs="i")

strps <- seq(1, nlevels(res$pmid)+1, 2)
rect(-100, strps-0.5, 100, strps+0.5, col = "gray92", border = NA)
abline(v = 0, lty = 3, col = "black")

points(jitter(log2(cap_or(res$diff_odds_neu))), jitter(as.numeric(res$pmid), 0.7), pch = ifelse(res$diff_pval_neu <= 0.05, 19, 1), col = cell2col("neu"))
points(jitter(log2(cap_or(res$diff_odds_mono))), jitter(as.numeric(res$pmid), 0.7), pch = ifelse(res$diff_pval_mono <= 0.05, 19, 1), col = cell2col("mono"))
points(jitter(log2(cap_or(res$diff_odds_tcd4m))), jitter(as.numeric(res$pmid), 0.7), pch = ifelse(res$diff_pval_tcd4m <= 0.05, 19, 1), col = cell2col("tcd4m"))
points(jitter(log2(cap_or(res$diff_odds_tcd4n))), jitter(as.numeric(res$pmid), 0.7), pch = ifelse(res$diff_pval_tcd4n <= 0.05, 19, 1), col = cell2col("tcd4n"))

axis(1, at = seq(-10,10,5))
axis(1, at = 13, labels = c("> 10"))
mtext("Log2 odds ratio", 1, at = 0, line = 1.5, cex = 0.6)


par(mar = c(3,0,2,1), tck = -0.02, mgp = c(2,0.5,0))
plot.new()
plot.window(xlim = c(0.5, 3.5), ylim = c(nlevels(res$pmid)+0.5, 0.5), xaxs = "i", yaxs = "i")

strps <- seq(1, nlevels(res$pmid)+1, 2)
rect(-100, strps-0.5, 100, strps+0.5, col = "gray92", border = NA)

text(1:3, 0.25, c("Odds ratio", "Predicted", "p-value"), pos = 3, xpd = TRUE, font = 2)
odds <- mapply(format_odds, res$diff_odds_tcd4n[!duplicated(res$pmid)])
text(1, 1:nlevels(res$pmid), ifelse(odds > 9999, expression(infinity), odds), cex = ifelse(odds > 9999, 1.75, 1))
text(2, 1:nlevels(res$pmid), mapply(format_perc, tapply(res$diff_perc_tcd4n, res$pmid, max)))
text(3, 1:nlevels(res$pmid), mapply(format_pvals, tapply(res$diff_pval_tcd4n, res$pmid, min)))

invisible(dev.off())
