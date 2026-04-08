#--- init ----------------------------------------------------------------------

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


#--- aging probes --------------------------------------------------------------

model  <- model.matrix(~ age + celltype, data = B$"")
model0 <- model.matrix(~ celltype, data = B$"")
res    <- row_lm_f(as.matrix(B), model, model0)

resB["pvalue.age"] <- res$pvalue


model  <- model.matrix(~ age + celltype, data = R$"")
model0 <- model.matrix(~ celltype, data = R$"")
res    <- row_lm_f(as.matrix(R), model, model0)

resR["pvalue.age"] <- res$pvalue


#--- celltype pvalues ----------------------------------------------------------

model  <- model.matrix(~ age + celltype, data = B$"")
model0 <- model.matrix(~ age, data = B$"")
res    <- row_lm_f(as.matrix(B), model, model0)

resB["pvalue.celltype"] <- res$pvalue


model  <- model.matrix(~ age + celltype, data = R$"")
model0 <- model.matrix(~ age, data = R$"")
res    <- row_lm_f(as.matrix(R), model, model0)

resR["pvalue.celltype"] <- res$pvalue


#--- stratify cytosines --------------------------------------------------------

is_diff  <- p.adjust(resB$pvalue, "fdr") <= 0.05
is_agect <- p.adjust(resB$pvalue.age, "fdr") <= 0.05 & p.adjust(resB$pvalue.celltype, "fdr") <= 0.05
is_age   <- p.adjust(resB$pvalue.age, "fdr") <= 0.05
is_ct    <- p.adjust(resB$pvalue.cell, "fdr") <= 0.05
resB$category <- ifelse(is_diff, "cell-dea",
                        ifelse(is_agect, "age-ct",
                               ifelse(is_age, "age-only",
                                      ifelse(is_ct, "ct-only", "rest"))))

is_diff  <- p.adjust(resR$pvalue, "fdr") <= 0.05
is_agect <- p.adjust(resR$pvalue.age, "fdr") <= 0.05 & p.adjust(resR$pvalue.celltype, "fdr") <= 0.05
is_age   <- p.adjust(resR$pvalue.age, "fdr") <= 0.05
is_ct    <- p.adjust(resR$pvalue.cell, "fdr") <= 0.05
resR$category <- ifelse(is_diff, "cell-dea",
                        ifelse(is_agect, "age-ct",
                               ifelse(is_age, "age-only",
                                      ifelse(is_ct, "ct-only", "rest"))))


#--- combine -------------------------------------------------------------------

cgs <- intersect(rownames(resB), rownames(resR))

resB <- resB[cgs,]
resR <- resR[cgs,]


res <- data.frame(id = cgs, row.names = cgs)

res$category <- ifelse(resB$category == "cell-dea" | resR$category == "cell-dea", "cell-dea",
                       ifelse(resB$category == "age-ct" | resR$category == "age-ct", "age-ct",
                              ifelse(resB$category == "age-only" | resR$category == "age-only", "age-only",
                                     ifelse(resB$category == "ct-only" | resR$category == "ct-only", "ct-only", "rest")
                       )))


#--- save ----------------------------------------------------------------------

saveRDS(res, "diffage_cgs_5categories.rds")

