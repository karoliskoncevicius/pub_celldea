#--- init ----------------------------------------------------------------------

library(annmatrix)
library(kklibrary)
library(matrixTests)

X <- readRDS("hannon/3_annmatrix.rds")


#--- diff ageing probes --------------------------------------------------------

model  <- model.matrix(~ -1 + age * celltype - age, data = X$"")
model0 <- model.matrix(~ -1 + age + celltype, data = X$"")
resH   <- row_lm_f(as.matrix(X), model, model0)

betanames <- sub("celltype", "", colnames(model))
betanames <- sub(":", ".", betanames)
colnames(resH)[grep("beta\\.", colnames(resH))] <- paste0("beta.", betanames)


#--- celltype pvalues ----------------------------------------------------------

resH["pvalue.celltype"] <- row_oneway_equalvar(X, X$celltype)$pvalue


#--- aging pvalues -------------------------------------------------------------

for(cl in unique(X$celltype)) {
  resH[paste0("pvalue.age.", cl)] <- row_cor_pearson(X[,X$celltype == cl], X$age[X$celltype == cl])$pvalue
}


#--- stratify cytosines --------------------------------------------------------

is_diff <- p.adjust(resH$pvalue, "fdr") <= 0.05
is_stat <- p.adjust(resH$pvalue.cell, "fdr") <= 0.05
resH$category <- ifelse(is_diff, "cell-dea", ifelse(is_stat, "age-static", "rest"))


#--- final ---------------------------------------------------------------------

res <- data.frame(id = rownames(resH), row.names = rownames(resH))

res$agingrate_neunpos   <- resH$beta.age.neunpos * 100
res$agingrate_sox10pos  <- resH$beta.age.sox10pos * 100
res$agingrate_doubleneg <- resH$beta.age.doubleneg * 100
res$pvalue_neunpos      <- resH$pvalue.age.neunpos
res$pvalue_sox10pos     <- resH$pvalue.age.sox10pos
res$pvalue_doubleneg    <- resH$pvalue.age.doubleneg
res$anova_pvalue        <- resH$pvalue.celltype
res$anova_qvalue        <- p.adjust(resH$pvalue.celltype, "fdr")
res$celldea_pvalue      <- resH$pvalue
res$celldea_qvalue      <- p.adjust(resH$pvalue, "fdr")
res$category            <- resH$category


#--- save ----------------------------------------------------------------------

saveRDS(res, "diffage_cgs_hannon.rds")

