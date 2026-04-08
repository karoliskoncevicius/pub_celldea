#--- init ----------------------------------------------------------------------

library(kklibrary)
library(annmatrix)
library(data.table)

key <- readRDS("samplekey.rds")


#--- read ----------------------------------------------------------------------

dat1 <- as.matrix(fread("input/GSE80417_normalizedBetas.csv"),  rownames = 1)
dat2 <- as.matrix(fread("input/GSE84727_normalisedBetas.csv"),  rownames = 1)
dat3 <- as.matrix(fread("input/GSE147221_Dublin_blood_processed_signals.csv"), rownames = 1)
dat4 <- as.matrix(fread("input/GSE152026_EUGEI_processed_signals.csv"), rownames = 1)
dat5 <- as.matrix(fread("input/GSE152027_IOP_processed_signals.csv"), rownames = 1)

# clean
dat3 <- dat3[,!grepl("_Detection", colnames(dat3))]
dat4 <- dat4[,!grepl("_Detection", colnames(dat4))]
dat5 <- dat5[,!grepl("_Detection", colnames(dat5))]

cpgs <- rownames(dat1) |> intersect(rownames(dat2)) |> intersect(rownames(dat3)) |> intersect(rownames(dat4)) |> intersect(rownames(dat5))
dat1 <- dat1[cpgs,]
dat2 <- dat2[cpgs,]
dat3 <- dat3[cpgs,]
dat4 <- dat4[cpgs,]
dat5 <- dat5[cpgs,]

dat <- cbind(dat1, dat2, dat3, dat4, dat5)
rm(dat1, dat2, dat3, dat4, dat5)

key <- key[colnames(dat),]
X   <- annmatrix(dat, cann = key)


#--- save ----------------------------------------------------------------------

saveRDS(X, "1_normalized.rds")
