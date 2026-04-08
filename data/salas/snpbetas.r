#--- init ----------------------------------------------------------------------

library(kklibrary)
suppressPackageStartupMessages(library(minfi))


#--- dependencies --------------------------------------------------------------

outfile <- commandArgs(TRUE)[3]

dat <- readRDS("1_raw.rds")


#--- get snp probes ------------------------------------------------------------

snps <- getSnpBeta(dat)


#--- save ----------------------------------------------------------------------

saveRDS(snps, "snpbetas.rds")

