#--- init ----------------------------------------------------------------------

library(kklibrary)

df1 <- read.soft("input/GSE80417.txt")
df2 <- read.soft("input/GSE84727.txt")
df3 <- read.soft("input/GSE147221.txt")
df4 <- read.soft("input/GSE152026.txt")
df5 <- read.soft("input/GSE152027.txt")


#--- prepare -------------------------------------------------------------------

key1 <- data.frame(id = df1$description,
                   cohort = "gse80417",
                   sentrix_id = sub("_.*", "", df1$description),
                   diagnosis = ifelse(df1$sample_disease_status == 1, "control", "scz"),
                   sex = df1$sample_Sex,
                   age = df1$sample_age
                   )

key2 <- data.frame(id = df2$sample_sentrixids,
                   cohort = "gse84727",
                   sentrix_id = sub("_.*", "", df2$sample_sentrixids),
                   diagnosis = ifelse(df2$sample_disease_status == 1, "control", "scz"),
                   sex = df2$sample_Sex,
                   age = df2$sample_age
                   )

key3 <- data.frame(id = sub(":.*", "", df3$title),
                   cohort = "gse147221",
                   sentrix_id = sub("_.*", "", sub(":.*", "", df3$title)),
                   diagnosis = ifelse(tolower(df3$sample_status) == "case", "scz", "control"),
                   sex = df3$sample_Sex,
                   age = df3$sample_age
                   )

key4 <- data.frame(id = sub(" .*", "", df4$title),
                   cohort = "gse152026",
                   sentrix_id = sub("_.*", "", sub(" .*", "", df4$title)),
                   diagnosis = ifelse(tolower(df4$sample_phenotype) == "case", "fep", "control"),
                   sex = df4$sample_Sex,
                   age = df4$sample_age
                   )

key5 <- data.frame(id = sub(" .*", "", df5$title),
                   cohort = "gse152027",
                   sentrix_id = sub("_.*", "", sub(" .*", "", df5$title)),
                   diagnosis = c(CON = "control", FEP = "fep", SCZ = "scz")[df5$sample_status],
                   sex = df5$sample_gender,
                   age = df5$sample_ageatbloodcollection
                   )


key <- rbind(key1, key2, key3, key4, key5)
rownames(key) <- key$id


#--- save ----------------------------------------------------------------------

saveRDS(key, "samplekey.rds")

