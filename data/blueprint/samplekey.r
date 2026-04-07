#--- init ----------------------------------------------------------------------

library(kklibrary)

df <- read.delim("input/samplekey.txt")


#--- prepare -------------------------------------------------------------------

key <- data.frame(id = df$ID)

key$sentrix_id   <- as.character(df$Sentrix_ID)
key$sentrix_row  <- substr(df$Sentrix_Position, 2, 3)
key$sentrix_col  <- substr(df$Sentrix_Position, 5, 6)

key$lab          <- df$Centre
key$batch        <- as.character(df$Sample_Batch)
key$well_row     <- sub("[0-9]+", "", df$Sample_Well)
key$well_col     <- sub("[A-Z]+", "", df$Sample_Well)

key$comment      <- df$Comment
key$is_repeated  <- df$Repeat == 1
key$is_crossover <- df$Cross_Over == 1

key$donor        <- df$Donor
key$celltype     <- c(NEU = "neu", MO = "mono", Tcell = "tcd4n")[df$Sample_Group]
key$sex          <- df$Sex
key$age          <- df$Age_at_visit

rownames(key) <- key$id

key <- key[order(key$id),]


#--- save ----------------------------------------------------------------------

saveRDS(key, "samplekey.rds")

