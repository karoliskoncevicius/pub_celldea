#--- init ----------------------------------------------------------------------

library(kklibrary)


#--- dependencies --------------------------------------------------------------

df <- read.csv("input/samplekey.csv", skip = 7)

# NOTE: fix one age (seems like a mistake for one cell type)
df$Age[which(df$Donor_ID == 2007)] <- 53


#--- prepare -------------------------------------------------------------------

key <- data.frame(id = paste0(df$Sentrix_ID, "_", df$Sentrix_Position))

key$sentrix_id  <- as.character(df$Sentrix_ID)
key$sentrix_row <- sub("R(..)C(..).*", "\\1", df$Sentrix_Position)
key$sentrix_col <- sub("R(..)C(..).*", "\\2", df$Sentrix_Position)

key$batch       <- df$Sample_Batch
key$pool_id     <- df$Pool_ID
key$well_row    <- substr(df$Sample_Well, 1, 1)
key$well_col    <- substr(df$Sample_Well, 2, 3)
key$facs_purity <- df$Cell_Purity_FACS

key$celltype[df$Sample_Group == "Monocyte"] <- "mono"
key$celltype[df$Sample_Group == "B_cell"]   <- "bcell"
key$celltype[df$Sample_Group == "T_cell"]   <- "tcd4"

key$diagnosis     <- ifelse(df$Status == "Y", "T1D", "control")
key$diagnosis_age <- df$Age_Diagnosis
key$bce           <- df$BCE
key$gad           <- df$GAD
key$ia2           <- df$IA2
key$znt8          <- df$ZnT8
key$donor         <- paste0("S", df$Donor_ID)
key$residence     <- df$Residence
key$sex           <- df$Sex
key$age           <- df$Age_Collection

rownames(key) <- key$id

key <- key[order(rownames(key)),]


#--- save ----------------------------------------------------------------------

saveRDS(key, "samplekey.rds")

