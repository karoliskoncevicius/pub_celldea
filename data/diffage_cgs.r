#!/usr/bin/env Rscript

library(annmatrix)
library(kklibrary)
library(matrixTests)

B <- readRDS("blueprint/3_annmatrix.rds")
R <- readRDS("rhead/3_annmatrix.rds")


#--- diff ageing probes --------------------------------------------------------

model  <- model.matrix(~ -1 + age * celltype - age, data = B$"")
model0 <- model.matrix(~ -1 + age + celltype, data = B$"")
resB   <- row_lm_f(as.matrix(B), model, model0)

betanames <- sub("celltype", "", colnames(model))
betanames <- sub(":", ".", betanames)
colnames(resB)[grep("beta\\.", colnames(resB))] <- paste0("beta.", betanames)


model  <- model.matrix(~ -1 + age * celltype - age, data = R$"")
model0 <- model.matrix(~ -1 + age + celltype, data = R$"")
resR   <- row_lm_f(as.matrix(R), model, model0)

betanames <- sub("celltype", "", colnames(model))
betanames <- sub(":", ".", betanames)
colnames(resR)[grep("beta\\.", colnames(resR))] <- paste0("beta.", betanames)


#--- celltype pvalues ----------------------------------------------------------

resB["pvalue.celltype"] <- row_oneway_equalvar(B, B$celltype)$pvalue

resR["pvalue.celltype"] <- row_oneway_equalvar(R, R$celltype)$pvalue


#--- aging pvalues -------------------------------------------------------------

for(cl in unique(B$celltype)) {
  resB[paste0("pvalue.age.", cl)] <- row_cor_pearson(B[,B$celltype == cl], B$age[B$celltype == cl])$pvalue
}

for(cl in unique(R$celltype)) {
  resR[paste0("pvalue.age.", cl)] <- row_cor_pearson(R[,R$celltype == cl], R$age[R$celltype == cl])$pvalue
}


#--- stratify cytosines --------------------------------------------------------

is_diff <- p.adjust(resB$pvalue, "fdr") <= 0.05
is_stat <- p.adjust(resB$pvalue.cell, "fdr") <= 0.05
resB$category <- ifelse(is_diff, "cell-dea", ifelse(is_stat, "age-static", "rest"))

is_diff <- p.adjust(resR$pvalue, "fdr") <= 0.05
is_stat <- p.adjust(resR$pvalue.cell, "fdr") <= 0.05
resR$category <- ifelse(is_diff, "cell-dea", ifelse(is_stat, "age-static", "rest"))


#--- combine -------------------------------------------------------------------

cgs <- intersect(rownames(resB), rownames(resR))

resB <- resB[cgs,]
resR <- resR[cgs,]


res <- data.frame(id = cgs, row.names = cgs)

res$bp_agingrate_neu   <- resB$beta.age.neu * 100
res$bp_agingrate_mono  <- resB$beta.age.mono * 100
res$bp_agingrate_tcd4n <- resB$beta.age.tcd4n * 100
res$bp_anova_pvalue    <- resB$pvalue.celltype
res$bp_anova_qvalue    <- p.adjust(resB$pvalue.celltype, "fdr")
res$bp_celldea_pvalue  <- resB$pvalue
res$bp_celldea_qvalue  <- p.adjust(resB$pvalue, "fdr")

res$ra_agingrate_mono  <- resR$beta.age.mono * 100
res$ra_agingrate_tcd4n <- resR$beta.age.tcd4n * 100
res$ra_agingrate_tcd4m <- resR$beta.age.tcd4m * 100
res$ra_anova_pvalue    <- resR$pvalue.celltype
res$ra_anova_qvalue    <- p.adjust(resR$pvalue.celltype, "fdr")
res$ra_celldea_pvalue  <- resR$pvalue
res$ra_celldea_qvalue  <- p.adjust(resR$pvalue, "fdr")

res$category <- ifelse(resB$category == "cell-dea" | resR$category == "cell-dea", "cell-dea",
                       ifelse(resB$category == "age-static" | resR$category == "age-static", "age-static", "rest")
                       )


#--- save ----------------------------------------------------------------------

saveRDS(res, "diffage_cgs.rds")

