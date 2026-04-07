#--- init ----------------------------------------------------------------------

library(kklibrary)
library(annmatrix)
suppressPackageStartupMessages(library(inops))
suppressPackageStartupMessages(library(minfi))

masks <- read.delim("input/EPIC.hg38.mask.tsv")

dat <- readRDS("2_normalized.rds")
key <- readRDS("samplekey.rds")
S   <- readRDS("snpbetas.rds")
D   <- readRDS("detectionpvalues.rds")

stopifnot(identical(colnames(dat), rownames(key)))
stopifnot(identical(colnames(dat), colnames(S)))
stopifnot(identical(colnames(dat), colnames(D)))


#--- construct annmatrix -------------------------------------------------------

X <- annmatrix(getBeta(dat), getAnnotation(dat), key)

names(X@'') <- tolower(names(X@''))
X@'' <- data.frame(id = rownames(X), X@'')


#--- filter cpgs ---------------------------------------------------------------

# ch probes
X <- X[!grepl("^ch|^rs|^nv", X@id), ]

# sex chromosomes
X <- X[X@chr %out% c("chrX", "chrY"), ]

# masked probes
X <- X[!(X@id %in% masks$probeID[masks$MASK_general]), ]


#--- remove outliers -----------------------------------------------------------

# samples older than 100
X <- X[,X$age <= 100]

# mixed brain samples and celltypes with low sample size
X <- X[,X$celltype %out% c("brain", "irf8pos", "tripleneg")]

# missing age info
X <- X[,!is.na(X$age)]

# bad sex info
dat <- dat[,colnames(X)]
X   <- X[,is.na(X$sex) | getSex(dat)$predictedSex == X$sex]

# poor detection
D <- D[rownames(X), colnames(X)]
X <- X[,colMeans(D > 0.01) < 0.05]

# inter array correlation
X <- X[,iac(X, X$celltype)[3,] > -3]


#--- mask values ---------------------------------------------------------------

# remove poor detection probes
D <- D[rownames(X), colnames(X)]
X <- X[rowMeans(D > 0.01) < 0.05,]

# mask remaining poor detection values
D <- D[rownames(X), colnames(X)]
X[D > 0.01] <- NA


#--- save ----------------------------------------------------------------------

saveRDS(X, "3_annmatrix.rds")

