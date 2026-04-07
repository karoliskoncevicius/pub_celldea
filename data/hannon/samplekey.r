#--- init ----------------------------------------------------------------------

library(kklibrary)


#--- dependencies --------------------------------------------------------------

df <- read.soft("input/GSE279509.txt")


#--- prepare -------------------------------------------------------------------

key <- data.frame(id = sub("_Grn.*", "", basename(df$supplementary_file)))

key$sentrix_id  <- sub(".*_(.*)_R(..)C(..).*", "\\1", df$supplementary_file)
key$sentrix_row <- sub(".*_(.*)_R(..)C(..).*", "\\2", df$supplementary_file)
key$sentrix_col <- sub(".*_(.*)_R(..)C(..).*", "\\3", df$supplementary_file)

key$celltype  <- tolower(df$sample_cell_type)
key$celltype[key$celltype == "total"] <- "brain"
key$antibody  <- df$sample_antibody_combination
key$diagnosis <- NA
key$donor     <- df$sample_individual_id
key$sex       <- df$sample_Sex
key$age       <- df$sample_age

rownames(key) <- key$id


#--- save ----------------------------------------------------------------------

saveRDS(key, "samplekey.rds")

