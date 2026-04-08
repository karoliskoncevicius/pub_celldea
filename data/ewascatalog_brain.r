#--- init ----------------------------------------------------------------------

suppressPackageStartupMessages(library(inops))

res <- read.delim("EWAScatalog/ewascatalog-results.txt")
std <- read.delim("EWAScatalog/ewascatalog-studies.txt")


#--- clean ---------------------------------------------------------------------

res <- res[,c("StudyID", "CpG", "P", "Beta")]

std$Tissue     <- tolower(std$Tissue)
std$Trait      <- tolower(std$Trait)
std$Covariates <- tolower(std$Covariates)
std$Age        <- tolower(std$Age)

std <- std[,c("StudyID", "Author", "PMID", "Date", "Trait", "Covariates", "N", "Methylation_Array", "Tissue", "Age")]


#--- filter --------------------------------------------------------------------

std <- std[std$Methylation_Array %out~% "sequencing",]
std <- std[std$Tissue %in~% c("brain", "temporal gyrus", "cortex", "cerebellum"),]
std <- std[!is.na(std$PMID) & std$PMID %in#% 1:10,]
std <- std[std$Trait %out~% c("^age", "^sex$"),]

std$Age               <- NULL
std$Tissue            <- NULL
std$Methylation_Array <- NULL

std$Year <- std$Date
std$Year[grepl("-", std$Year)]      <- sub("-.*", "", std$Year[grepl("-", std$Year)])
std$Year[grepl("^[A-Z]", std$Year)] <- sub(".*, ", "", std$Year[grepl("^[A-Z]", std$Year)])

std$Author <- sub("^De ", "", std$Author)
std$Author <- sub(" .*", "", std$Author)


#--- merge ---------------------------------------------------------------------

res  <- res[res$StudyID %in% std$StudyID, ]
ewas <- data.frame(res, std[match(res$StudyID, std$StudyID),-1] )

ewas <- ewas[grepl("cell composition", ewas$Covariates),]
ewas <- ewas[!grepl("cell composition \\(reference free\\)", ewas$Covariates),]


#--- reorder -------------------------------------------------------------------

ewas <- ewas[,c("StudyID", "PMID", "Author", "Year", "Trait", "Covariates", "N", "CpG", "P", "Beta")]
colnames(ewas) <- tolower(colnames(ewas))


#--- save ----------------------------------------------------------------------

saveRDS(ewas, "ewascatalog_brain.rds")

