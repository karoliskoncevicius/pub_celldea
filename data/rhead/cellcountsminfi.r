#--- init ----------------------------------------------------------------------

library(kklibrary)
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(FlowSorted.Blood.EPIC))
suppressPackageStartupMessages(library(FlowSorted.BloodExtended.EPIC))

X <- readRDS("1_raw.rds")


#--- estimate cell types -------------------------------------------------------

cells <- estimateCellCounts2(X, compositeCellType = "BloodExtended", processMethod = "preprocessNoob", probeSelect = "both",
                             cellTypes = c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")
                             )$prop

# rename
newnames <- c(Neu = "neu", Mono = "mono", Eos = "eos", Bas = "baso", NK = "nk", Bnv = "bn", Bmem = "bm", Treg = "treg", CD4nv = "tcd4n", CD4mem = "tcd4m", CD8nv = "tcd8n", CD8mem = "tcd8m")

colnames(cells) <- newnames[colnames(cells)]  # rename
cells           <- cells[,newnames]           # reorder


#--- save ----------------------------------------------------------------------

saveRDS(cells, "cellcountsminfi.rds")

