#--- init ----------------------------------------------------------------------

library(EpiDISH)
library(basetheme)
library(kklibrary)
library(annmatrix)

B   <- readRDS("data/blueprint/3_annmatrix.rds")
R   <- readRDS("data/rhead/3_annmatrix.rds")
dfa <- readRDS("data/diffage_cgs.rds")

ref      <- cent12CT450k.m
dea_cpgs <- rownames(dfa)[dfa$category == "cell-dea"]

# CD4Tnv=1, CD4Tmem=3, CD8Tmem=7, CD8Tnv=8, Treg=6
cols <- c(1, 3, 6, 7, 8)
cpgs <- unlist(lapply(cols, function(b) rownames(ref)[(b-1)*50 + 1:50]))

dea_in_t    <- cpgs[cpgs %in% dea_cpgs]
nondea_in_t <- cpgs[!cpgs %in% dea_cpgs]

set.seed(1)
dea_in_t    <- sample(dea_in_t)
nondea_in_t <- sample(nondea_in_t)

run_analysis <- function(X, remove_cpgs) {
  ref_iter <- ref
  pool     <- remove_cpgs[remove_cpgs %in% rownames(ref_iter)]
  n        <- length(pool)

  inds <- X$celltype == "tcd4n"
  ages <- X$age[inds]

  cors <- pvals <- numeric(n + 1)

  for(i in 0:n) {
    est        <- epidish(X, ref_iter, method="RPC")$estF
    ct         <- cor.test(ages, est[inds, "CD4Tnv"])
    cors[i+1]  <- ct$estimate
    pvals[i+1] <- ct$p.value
    if(i < n)
      ref_iter <- ref_iter[rownames(ref_iter) != pool[i+1], ]
  }
  list(n_removed=0:n, cor=cors, pval=pvals)
}

# run
results <- list()
results$B_dea    <- run_analysis(B, dea_in_t)
results$R_dea    <- run_analysis(R, dea_in_t)
results$B_nondea <- run_analysis(B, nondea_in_t)
results$R_nondea <- run_analysis(R, nondea_in_t)

# plot
pdf("fig_s9.pdf", width = 14/2.54, height = 4/2.54, pointsize = 8)
par(mfrow = c(1,4), mar = c(4,3,3,0), oma = c(0,0,0,1), tck = -0.02, mgp = c(2,0.5,0), las = 1)

ylim <- range(sapply(results, function(x) x$cor), na.rm=TRUE)

panels <- list(
  list(ds="B",     type="dea",    letter="a", title="Blueprint",    group_title="Removing cell-DEA cytosines"),
  list(ds="R", type="dea",    letter="",  title="RA",    group_title=NULL),
  list(ds="B",     type="nondea", letter="b", title="Blueprint",    group_title="Removing non-cell-DEA cytosines"),
  list(ds="R", type="nondea", letter="",  title="RA",    group_title=NULL)
)

for(p in panels) {
  res <- results[[paste(p$ds, p$type, sep="_")]]
  plot(res$n_removed,
       res$cor,
       pch=19,
       cex=0.8,
       col=cell2col("tcd4n"),
       xlim=c(0, max(res$n_removed)+2),
       ylim=ylim,
       xlab="",
       ylab="Correlation with age", main="", las=1, xaxt="n")

  title(xlab=if(p$type=="dea") "Number of\ncell-dea cytosines removed" else "Number of\nnon-cell-dea cytosines removed", line = 2.5)
  axis(1, at=seq(0,150,25))
  title(p$title, line=0.5, adj=0)
  # if(!is.null(p$group_title))
    # title(p$group_title, line = 1.5, adj = 0, xpd = NA)
  if(p$letter != "")
    mtext(p$letter, 3, adj= -0.15, line=1.5, font=2)
}

invisible(dev.off())
