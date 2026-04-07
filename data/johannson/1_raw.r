#--- init ----------------------------------------------------------------------

suppressPackageStartupMessages(library(minfi))


#--- dependencies --------------------------------------------------------------

infiles <- list.files("input", pattern = "Grn.idat$", recursive = TRUE, full.names = TRUE)


#--- read ----------------------------------------------------------------------

dat <- read.metharray(infiles, force = TRUE)


#--- save ----------------------------------------------------------------------

saveRDS(dat, "1_raw.rds")

