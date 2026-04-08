#--- init ----------------------------------------------------------------------

library(annmatrix)
library(kklibrary)
library(matrixTests)

X <- readRDS("rhead/3_annmatrix.rds")


#--- diff ageing probes --------------------------------------------------------

model  <- model.matrix(~ -1 + age * celltype - age, data = X$"")
model0 <- model.matrix(~ -1 + age + celltype, data = X$"")
resB   <- row_lm_f(as.matrix(X), model, model0)

betanames <- sub("celltype", "", colnames(model))
betanames <- sub(":", ".", betanames)
colnames(resB)[grep("beta\\.", colnames(resB))] <- paste0("beta.", betanames)


#--- celltype pvalues ----------------------------------------------------------

resB["pvalue.celltype"] <- row_oneway_equalvar(X, X$celltype)$pvalue


#--- aging pvalues -------------------------------------------------------------

for(cl in unique(X$celltype)) {
  resB[paste0("pvalue.age.", cl)] <- row_cor_pearson(X[,X$celltype == cl], X$age[X$celltype == cl])$pvalue
}


#--- stratify cytosines --------------------------------------------------------

is_diff <- p.adjust(resB$pvalue, "fdr") <= 0.05
is_stat <- p.adjust(resB$pvalue.cell, "fdr") <= 0.05
resB$category <- ifelse(is_diff, "cell-dea", ifelse(is_stat, "age-static", "rest"))


#--- final ---------------------------------------------------------------------

res <- data.frame(id = rownames(resB), row.names = rownames(resB))

res$agingrate_mono  <- resB$beta.age.mono * 100
res$agingrate_tcd4n <- resB$beta.age.tcd4n * 100
res$agingrate_tcd4m <- resB$beta.age.tcd4m * 100
res$pvalue_mono     <- resB$pvalue.age.mono
res$pvalue_tcd4n    <- resB$pvalue.age.tcd4n
res$pvalue_tcd4m    <- resB$pvalue.age.tcd4m
res$anova_pvalue    <- resB$pvalue.celltype
res$anova_qvalue    <- p.adjust(resB$pvalue.celltype, "fdr")
res$celldea_pvalue  <- resB$pvalue
res$celldea_qvalue  <- p.adjust(resB$pvalue, "fdr")
res$category        <- resB$category


#--- save ----------------------------------------------------------------------

saveRDS(res, "diffage_cgs_rhead.rds")

