#--- init ----------------------------------------------------------------------

library(kklibrary)
suppressPackageStartupMessages(library(minfi))


dat <- readRDS("1_raw.rds")


#--- normalize -----------------------------------------------------------------

dat <- mapToGenome(preprocessNoob(dat))


#--- save ----------------------------------------------------------------------

saveRDS(dat, "2_normalized.rds")

