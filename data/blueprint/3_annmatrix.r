#--- init ----------------------------------------------------------------------

library(kklibrary)
library(annmatrix)
suppressPackageStartupMessages(library(inops))
suppressPackageStartupMessages(library(minfi))


#--- dependencies --------------------------------------------------------------

masks <- read.delim("input/HM450.hg38.mask.tsv")

dat <- readRDS("2_normalized.rds")
key <- readRDS("samplekey.rds")
S   <- readRDS("snpbetas.rds")
D   <- readRDS("detectionpvalues.rds")
C   <- readRDS("cellcountsepidish.rds")

stopifnot(identical(colnames(dat), rownames(key)))
stopifnot(identical(colnames(dat), colnames(S)))
stopifnot(identical(colnames(dat), colnames(D)))
stopifnot(identical(colnames(dat), rownames(C)))


#--- construct annmatrix -------------------------------------------------------

X <- annmatrix(getBeta(dat), getAnnotation(dat), key)

names(X@'') <- tolower(names(X@''))
X@'' <- data.frame(id = rownames(X), X@'')


#--- filter cpgs ---------------------------------------------------------------

# non-cg probes
X <- X[X@id %out~% c("^ch|^rs|^nv"), ]

# sex chromosomes
X <- X[X@chr %out% c("chrX", "chrY"), ]

# masked probes
X <- X[X@id %out% masks$probeID[masks$MASK_general], ]


#--- remove outliers -----------------------------------------------------------

# known issues
X <- X[,is.na(X$comment)]
X$comment <- NULL

# repeat or cross-over
X <- X[,!X$is_repeated & !X$is_crossover]
X$is_repeated <- X$is_crossover <- NULL

# bad sex info
dat <- dat[,colnames(X)]
X   <- X[,is.na(X$sex) | getSex(dat)$predictedSex == X$sex]

# poor detection
D <- D[rownames(X), colnames(X)]
X <- X[,colMeans(D > 0.01) < 0.01]

# poor purity
C   <- C[colnames(X),]
pur <- C[cbind(colnames(X), X$celltype)]
X   <- X[,ave(pur, X$celltype, FUN = scale) > -3]

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

