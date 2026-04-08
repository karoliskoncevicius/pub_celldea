#--- init ----------------------------------------------------------------------

library(kklibrary)

df <- read.soft("input/GSE315327.txt")


#--- prepare -------------------------------------------------------------------

key <- data.frame(id = df$sample_sentrix_pos)

key$sentrix_id  <- sub(".*_(.*)_R(..)C(..).*", "\\1", df$supplementary_file)
key$sentrix_row <- sub(".*_(.*)_R(..)C(..).*", "\\2", df$supplementary_file)
key$sentrix_col <- sub(".*_(.*)_R(..)C(..).*", "\\3", df$supplementary_file)

key$well_row    <- sub("[0-9]+", "", df$sample_well_pos)
key$well_col    <- sub("[A-~]+", "", df$sample_well_pos)

key$celltype    <- ifelse(startsWith(df$sample_cell_type, "n"), "neu", "tcell")
key$donor       <- df$sample_donor

key$smoking     <- df$sample_smoking
key$bmi         <- df$sample_bmi
key$diagnosis   <- df$sample_diagnosis
key$medication  <- df$sample_medication
key$sex         <- df$sample_Sex
key$age         <- df$sample_age

key$purity      <- df$sample_purity
key$perc_tcd4   <- df$sample_perc_tcd4
key$perc_tcd8   <- df$sample_perc_tcd8
key$perc_tcd4n  <- df$sample_perc_tcd4n
key$perc_tcd4m  <- df$sample_perc_tcd4m

rownames(key) <- key$id
key <- key[order(key$id),]

#--- save ----------------------------------------------------------------------

saveRDS(key, "samplekey.rds")

