#--- init ----------------------------------------------------------------------

library(kklibrary)
library(EpiDISH)
suppressPackageStartupMessages(library(minfi))

X <- readRDS("2_normalized.rds")


#--- estimate cell types -------------------------------------------------------

# select reference
ref <- list(IlluminaHumanMethylation450k   = cent12CT450k.m,
            IlluminaHumanMethylationEPIC   = cent12CT.m,
            IlluminaHumanMethylationEPICv2 = cent12CT.m
            )[[annotation(X)["array"]]]

cells <- epidish(getBeta(X), ref)$estF

# rename
newnames <- c(Neu = "neu", Mono = "mono", Eos = "eos", Baso = "baso", NK = "nk", Bnv = "bn", Bmem = "bm", Treg = "treg", CD4Tnv = "tcd4n", CD4Tmem = "tcd4m", CD8Tnv = "tcd8n", CD8Tmem = "tcd8m")

colnames(cells) <- newnames[colnames(cells)]  # rename
cells           <- cells[,newnames]           # reorder


#--- save ----------------------------------------------------------------------

saveRDS(cells, "cellcountsepidish.rds")

