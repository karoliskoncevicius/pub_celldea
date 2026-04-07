#--- init ----------------------------------------------------------------------

library(kklibrary)
suppressPackageStartupMessages(library(minfi))


#--- dependencies --------------------------------------------------------------

dat <- readRDS("1_raw.rds")


#--- normalize -----------------------------------------------------------------

dat <- mapToGenome(preprocessNoob(dat))


#--- save ----------------------------------------------------------------------

saveRDS(dat, "2_normalized.rds")

