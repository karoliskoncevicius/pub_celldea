#--- init ----------------------------------------------------------------------

library(kklibrary)
library(ewastools)
suppressPackageStartupMessages(library(minfi))


#--- dependencies --------------------------------------------------------------

dat <- readRDS("1_raw.rds")


#--- get detection p-values ----------------------------------------------------

ps <- detectionP.minfi(dat)


#--- save ----------------------------------------------------------------------

saveRDS(ps, "detectionpvalues.rds")

