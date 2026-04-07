#--- init ----------------------------------------------------------------------

library(kklibrary)

df <- read.soft("input/GSE131989.txt")


#--- prepare -------------------------------------------------------------------

key <- data.frame(id = sub("_Grn.*", "", basename(df$supplementary_file)))

key$sentrix_id  <- sub(".*_(.*)_R(..)C(..).*", "\\1", df$supplementary_file)
key$sentrix_row <- sub(".*_(.*)_R(..)C(..).*", "\\2", df$supplementary_file)
key$sentrix_col <- sub(".*_(.*)_R(..)C(..).*", "\\3", df$supplementary_file)

key$rundate   <- df$sample_run_date
key$plate     <- df$sample_sample_plate
key$well_row  <- sub("[0-9]+", "", df$sample_sample_well)
key$well_col  <- sub("[A-Z]+", "", df$sample_sample_well)
key$issues    <- sub(".*\\((.*)\\)|.*", "\\1", df$title)

key$celltype  <- c(CD19 = "bcell", CD4Mem = "tcd4m", CD4Nv = "tcd4n", CD14 = "mono")[df$sample_sample_type]
key$race      <- tolower(df$sample_ethnicity)
key$diagnosis <- c("ra", "control")[df$sample_subjecttype]
key$smoking   <- as.logical(df$sample_smokerever)
key$sex       <- c(Female = "F", Male = "M")[df$sample_gender]
key$age       <- df$sample_ageatblooddraw

rownames(key) <- key$id


#--- save ----------------------------------------------------------------------

saveRDS(key, "samplekey.rds")

