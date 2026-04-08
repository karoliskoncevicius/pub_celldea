#--- init ----------------------------------------------------------------------

library(kklibrary)
library(ewastools)
suppressPackageStartupMessages(library(minfi))


#--- dependencies --------------------------------------------------------------

outfile <- commandArgs(TRUE)[3]

dat <- readRDS("1_raw.rds")


#--- get detection p-values ----------------------------------------------------

ps <- detectionP.minfi(dat)


#--- save ----------------------------------------------------------------------

saveRDS(ps, "detectionpvalues.rds")

