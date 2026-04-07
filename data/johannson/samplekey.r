#--- init ----------------------------------------------------------------------

library(kklibrary)

df <- read.soft("input/GSE87571.txt")


#--- prepare -------------------------------------------------------------------

key <- data.frame(id = sub("_Grn.*", "", basename(df$supplementary_file)))

key$sentrix_id  <- sub(".*_(.*)_R(..)C(..).*", "\\1", df$supplementary_file)
key$sentrix_row <- sub(".*_(.*)_R(..)C(..).*", "\\2", df$supplementary_file)
key$sentrix_col <- sub(".*_(.*)_R(..)C(..).*", "\\3", df$supplementary_file)

key$celltype  <- "wb"
key$diagnosis <- "control"
key$sex       <- substr(df$sample_gender, 1, 1)
key$age       <- df$sample_age

rownames(key) <- key$id


#--- save ----------------------------------------------------------------------

saveRDS(key, "samplekey.rds")

