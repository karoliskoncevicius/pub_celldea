#--- init ----------------------------------------------------------------------

library(readxl)
library(annmatrix)
library(kklibrary)
library(matrixTests)
library(EpiDISH)
library(DEAfree)
library(methylCIPHER)

X    <- readRDS("hannon_scz/1_normalized.rds")
df   <- readRDS("diffage_cgs.rds")

cent <- centDEAfree12CT.m


#--- clean ---------------------------------------------------------------------

# pre-select samples
X <- X[,!is.na(X$sex) & !is.na(X$age) & X$age < 100 & X$cohort %in% c("gse147221", "gse80417", "gse84727")]

# smoking score
X$smokingscore <- calcSmokingMcCartney(t(X))

# cell proportions
newnames <- c(Neu = "neu", Mono = "mono", Eos = "eos", Baso = "baso", NK = "nk", Bnv = "bn", Bmem = "bm", Treg = "treg", CD4Tnv = "tcd4n", CD4Tmem = "tcd4m", CD8Tnv = "tcd8n", CD8Tmem = "tcd8m")

centOld <- cent12CT450k.m
colnames(centOld) <- newnames[colnames(centOld)]

centOld <- centOld[,newnames] # reorder
cent    <- cent[,newnames]    # reorder

C  <- epidish(X, centDHSbloodDMC.m)$estF
C1 <- epidish(X, centOld)$estF
C2 <- epidish(X, cent)$estF

# intersect with diff-aging datasets
cgs <- intersect(rownames(df), rownames(X))
X   <- X[cgs,]


#--- results -------------------------------------------------------------------

# using 6 cell types
model  <- model.matrix(~ X$diagnosis + X$age + X$sex + X$cohort + X$smokingscore + C)
model0 <- model.matrix(~ X$age + X$sex + X$cohort + X$smokingscore + C)
res    <- row_lm_f(X, model, model0)

# adding interaction with age
model  <- model.matrix(~ X$diagnosis + X$age + X$sex + X$cohort + X$smokingscore + C*X$age)
model0 <- model.matrix(~ X$age + X$sex + X$cohort + X$smokingscore + C*X$age)
res2   <- row_lm_f(X, model, model0)

# using 12 cell types
model  <- model.matrix(~ X$diagnosis + X$age + X$sex + X$cohort + X$smokingscore + C1)
model0 <- model.matrix(~ X$age + X$sex + X$cohort + X$smokingscore + C1)
res3   <- row_lm_f(X, model, model0)

# using 12 cell types with age interaction
model  <- model.matrix(~ X$diagnosis + X$age + X$sex + X$cohort + X$smokingscore + C1*X$age)
model0 <- model.matrix(~ X$age + X$sex + X$cohort + X$smokingscore + C1*X$age)
res4   <- row_lm_f(X, model, model0)

# using new reference
model  <- model.matrix(~ X$diagnosis + X$age + X$sex + X$cohort + X$smokingscore + C2)
model0 <- model.matrix(~ X$age + X$sex + X$cohort + X$smokingscore + C2)
res5   <- row_lm_f(X, model, model0)

# using new reference with age interaction
model  <- model.matrix(~ X$diagnosis + X$age + X$sex + X$cohort + X$smokingscore + C2*X$age)
model0 <- model.matrix(~ X$age + X$sex + X$cohort + X$smokingscore + C2*X$age)
res6   <- row_lm_f(X, model, model0)


#--- save ----------------------------------------------------------------------

res <- list(res, res2, res3, res4, res5, res6)
names(res) <- c("epi_6ct", "epi_6ct_age", "epi_12ct", "epi_12ct_age", "new_12ct", "new_12ct_age")

saveRDS(res, "scz_ewas_correction.rds")
