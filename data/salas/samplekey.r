#--- init ----------------------------------------------------------------------

library(kklibrary)

df <- read.soft("input/GSE167998.txt")


#--- prepare -------------------------------------------------------------------

key <- data.frame(id = sub("_Grn.*", "", basename(df$supplementary_file)))

key$sentrix_id  <- sub(".*_(.*)_R(..)C(..).*", "\\1", df$supplementary_file)
key$sentrix_row <- sub(".*_(.*)_R(..)C(..).*", "\\2", df$supplementary_file)
key$sentrix_col <- sub(".*_(.*)_R(..)C(..).*", "\\3", df$supplementary_file)

key$purity     <- df$sample_purity
key$perc_bcell <- df$sample_bcell
key$perc_bn    <- df$sample_bnv
key$perc_bm    <- df$sample_bmem
key$perc_tcd4  <- df$sample_cd4t
key$perc_tcd4n <- df$sample_cd4nv
key$perc_tcd4m <- df$sample_cd4mem
key$perc_tcd8  <- df$sample_cd8t
key$perc_treg  <- df$sample_treg
key$perc_nk    <- df$sample_nk
key$perc_gran  <- df$sample_gran
key$perc_neu   <- df$sample_neu
key$perc_mono  <- df$sample_mono
key$perc_eos   <- df$sample_eos
key$issues     <- df$sample_flag_reason

key$celltype  <- c(Bas = "baso", Bmem = "bm", Bnv = "bn", CD4mem = "tcd4m", CD4nv = "tcd4n", CD8mem = "tcd8m", CD8nv = "tcd8n", Eos = "eos", MIX = "mix", Mono = "mono", Neu = "neu", NK = "nk", Treg = "treg")[df$description]
key$race      <- df$sample_ethnicity_wide
key$diagnosis <- "control"
key$smoking   <- df$sample_smoker == "Yes"
key$height    <- df$sample_height_cm
key$weight    <- df$sample_weight_kg
key$bmi       <- df$sample_bmi
key$age       <- df$sample_age

rownames(key) <- key$id


#--- save ----------------------------------------------------------------------

saveRDS(key, "samplekey.rds")

