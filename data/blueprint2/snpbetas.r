#--- init ----------------------------------------------------------------------

library(kklibrary)
suppressPackageStartupMessages(library(minfi))

dat <- readRDS("1_raw.rds")


#--- get snp probes ------------------------------------------------------------

snps <- getSnpBeta(dat)


#--- save ----------------------------------------------------------------------

saveRDS(snps, "snpbetas.rds")

