#!/usr/bin/env Rscript

library(annmatrix)
suppressPackageStartupMessages(library(inops))
suppressPackageStartupMessages(library(readxl))

res <- read.delim("input/ewascatalog-results.txt")
std <- read.delim("input/ewascatalog-studies.txt")


#--- clean ---------------------------------------------------------------------

res <- res[,c("StudyID", "CpG", "P", "Beta")]

std$Tissue     <- tolower(std$Tissue)
std$Trait      <- tolower(std$Trait)
std$Covariates <- tolower(std$Covariates)
std$Age        <- tolower(std$Age)

std <- std[,c("StudyID", "Author", "PMID", "Date", "Trait", "Covariates", "N", "Methylation_Array", "Tissue", "Age")]


#--- filter --------------------------------------------------------------------

std <- std[std$Methylation_Array %out~% "sequencing",]
std <- std[std$Tissue %in% c("blood", "leukocytes", "peripheral blood", "peripheral blood mononuclear cells", "whole blood"),]
std <- std[!is.na(std$PMID) & std$PMID %in#% 1:10,]
std <- std[std$Trait %out~% c("^age", "^aging", "^sex$"),]

std$PMID == "doi.org/10.1101/2023.11.02.23298000" <- "38692281"

std$Tissue            <- NULL
std$Methylation_Array <- NULL

std$Year <- std$Date
std$Year[grepl("-", std$Year)]      <- sub("-.*", "", std$Year[grepl("-", std$Year)])
std$Year[grepl("^[A-Z]", std$Year)] <- sub(".*, ", "", std$Year[grepl("^[A-Z]", std$Year)])

std$Author[std$Author == "Anna Ulrich"]        <- "Ulrich"
std$Author[std$Author == "Robert F. Hillary"]  <- "Hillary"
std$Author[std$Author == "Jiahui Si"]          <- "Si"
std$Author[std$Author == "Xueyi Shen"]         <- "Shen"
std$Author[std$Author == "Priyanka Choudhary"] <- "Choudhary"
std$Author[std$Author == "Jenny van Dongen"]   <- "Dongen"
std$Author[std$Author == "Ding, Yuan Chun"]    <- "Ding"
std$Author[std$Author == "Van Roekel E"]       <- "Roekel"
std$Author[std$Author == "Jiantao Ma"]         <- "Ma"
std$Author <- sub("^van ", "", std$Author)
std$Author <- sub("^der ", "", std$Author)
std$Author <- sub(" .*", "", std$Author)


#--- merge ---------------------------------------------------------------------

res  <- res[res$StudyID %in% std$StudyID, ]
ewas <- data.frame(res, std[match(res$StudyID, std$StudyID),-1] )

ewas <- ewas[grepl("cell composition", ewas$Covariates),]
ewas <- ewas[!grepl("cell composition \\(reference free\\)", ewas$Covariates),]


#--- reorder -------------------------------------------------------------------

ewas <- ewas[,c("StudyID", "PMID", "Author", "Year", "Trait", "Age", "Covariates", "N", "CpG", "P", "Beta")]
colnames(ewas) <- tolower(colnames(ewas))


#----- add ewas manually -------------------------------------------------------

# NOTE: adding additional EWAS that adjusted for 12 cell types (reportedly)

e <- as.data.frame(read_excel("input/KEPI_A_2289786_SM2258.xls", skip = 5))
e <- data.frame(studyid = "38090774_crohns", pmid = "38090774", author = "Reynolds", year = "2024", trait = "crohn's disease", age = "pediatric", covariates = NA, n = NA, cpg = e$CpGID, p = NA, beta = NA)
ewas <- rbind(ewas, e)

e <- as.data.frame(read_excel("input/42003_2026_9520_MOESM4_ESM.xlsx", sheet = 2, skip = 3))
e <- data.frame(studyid = "41639309_interleukin", pmid = "41639309", author = "Sinke", year = "2026", trait = "interleukin-6", age = "adult", covariates = NA, n = NA, cpg = e$ProbeID, p = e[,"padj (cells)"], beta = e[,"Beta (cells)"])
e <- e[e$p < 0.05,]
ewas <- rbind(ewas, e)


#--- save ----------------------------------------------------------------------

saveRDS(ewas, "ewascatalog_blood.rds")

