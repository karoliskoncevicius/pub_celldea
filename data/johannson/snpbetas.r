#--- inite ---------------------------------------------------------------------

library(kklibrary)
suppressPackageStartupMessages(library(minfi))


#--- dependencies --------------------------------------------------------------

dat <- readRDS("1_raw.rds")


#--- get snp probes ------------------------------------------------------------

snps <- getSnpBeta(dat)


#--- save ----------------------------------------------------------------------

saveRDS(snps, "snpbetas.rds")

