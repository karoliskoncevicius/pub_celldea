#--- init ----------------------------------------------------------------------

library(kklibrary)

dfage <- readRDS("data/diffage_cgs_blueprint.rds")


#--- age direction overlaps ----------------------------------------------------

dfage <- dfage[dfage$celldea_qvalue <= 0.05,]

sigs <- t(apply(dfage[,5:7], 1, function(x) p.adjust(x, "fdr"))) <= 0.05

agedir <- ifelse(sigs[,"pvalue_neu"], "neu", "")
agedir <- paste0(agedir, ifelse(sigs[,"pvalue_mono"], " mono", ""))
agedir <- paste0(agedir, ifelse(sigs[,"pvalue_tcd4n"], " tcd4n", ""))
agedir <- trimws(agedir)
agedir[agedir == ""] <- "none"

tab <- table(agedir)
tab <- sort(tab, decreasing = TRUE)


#--- plot ----------------------------------------------------------------------

pdf("fig_s1_a.pdf", width = 14/2.54, height = 6/2.54, pointsize = 8)

par(mar = c(1,5,1,1), tck=-0.01, mgp=c(2,0.5,0))

plot.new()
plot.window(xlim = c(0, 9), ylim = c(-11000,20000), xaxs = "i", yaxs = "i")

rect((1:8) - 0.46, 0, (1:8) + 0.45, tab, lwd = 2)
axis(2, at = seq(0,20000,5000), las = 1)
title(ylab = "Number of cytosines", line = 3.5, adj = 0.85)

cols <- c("white", cell2col("neu"))[grepl("neu", names(tab))+1]
cols <- c(cols, c("white", cell2col("mono"))[grepl("mono", names(tab))+1])
cols <- c(cols, c("white", cell2col("tcd4n"))[grepl("tcd4n", names(tab))+1])
symbols(rep(1:8, 3), rep(-c(3000,6000,9000), each = 8), circles = rep(0.19, 3*8), bg = cols, add = TRUE, inches = FALSE)
text(0, -c(3000,6000,9000), c("Neu", "Mono", "TCD4n"), pos = 2, xpd = NA)

invisible(dev.off())
