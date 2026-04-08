#--- init ----------------------------------------------------------------------

library(kklibrary)
suppressPackageStartupMessages(library(inops))

dfage <- readRDS("diffage_cgs_blueprint.rds")
ewas  <- readRDS("ewascatalog_blood.rds")


#--- filter --------------------------------------------------------------------

# only cpgs that are in dfage
ewas <- ewas[ewas$cpg %in% rownames(dfage),]

# only ewas with 50 or more cpgs
ewas <- ewas[ewas$studyid %out#% 1:49,]


#--- overlaps ------------------------------------------------------------------

is_stat <- dfage$category == "age-static"
is_diff <- dfage$category == "cell-dea"
is_rest <- dfage$category == "rest"

mytable <- function(x, y) {
  tab <- table(x, y)
  if (nrow(tab) == 1) {
    tab <- rbind(tab, c(0,0))
  }
  tab
}

res <- data.frame(study = unique(ewas$studyid))
for(s in unique(ewas$studyid)) {
  # static
  fish <- fisher.test(mytable(rownames(dfage) %in% ewas$cpg[ewas$studyid == s], is_stat))
  res$odds_cell[res$study == s] <- fish$estimate
  res$pval_cell[res$study == s] <- fish$p.value
  res$perc_cell[res$study == s] <- mean(is_stat[rownames(dfage) %in% ewas$cpg[ewas$studyid == s]]) * 100
  # diffage
  fish <- fisher.test(mytable(rownames(dfage) %in% ewas$cpg[ewas$studyid == s], is_diff))
  res$odds_diff[res$study == s] <- fish$estimate
  res$pval_diff[res$study == s] <- fish$p.value
  res$perc_diff[res$study == s] <- mean(is_diff[rownames(dfage) %in% ewas$cpg[ewas$studyid == s]]) * 100
  # rest
  fish <- fisher.test(mytable(rownames(dfage) %in% ewas$cpg[ewas$studyid == s], is_rest))
  res$odds_rest[res$study == s] <- fish$estimate
  res$pval_rest[res$study == s] <- fish$p.value
  res$perc_rest[res$study == s] <- mean(is_rest[rownames(dfage) %in% ewas$cpg[ewas$studyid == s]]) * 100
}


#--- add info ------------------------------------------------------------------

res$pmid   <- ewas$pmid[match(res$study, ewas$studyid)]
res$author <- ewas$author[match(res$study, ewas$studyid)]
res$year   <- ewas$year[match(res$study, ewas$studyid)]
res$trait  <- ewas$trait[match(res$study, ewas$studyid)]

res$trait %in~% "diet quality"           <- "diet quality"
res$trait %in~% "systemic lupus"         <- "systemic lupus erythematosus"
res$trait %in~% "sjogrens syndrome"      <- "sjogren's syndrome"
res$trait %in~% "cholesterol"            <- "serum cholesterol"
res$trait %in~% "cognitive abilities"    <- "cognitive abilities"
res$trait %in~% "alcohol consumption"    <- "alcohol consumption"
res$trait %in~% "glomerular"             <- "glomerular filtration rate"
res$trait %in~% "hypertension"           <- "pulmonary arterial hypertension"
res$trait %in~% "hepatitis"              <- "drug use and hepatitis"
res$trait %in~% "c-reactive protein"     <- "c-reactive protein"
res$trait %in~% "maternal smoking"       <- "maternal smoking"
res$trait %in~% "maternal education"     <- "maternal education"
res$trait %in% "bmi"                     <- "body mass index"
res$trait %in% "fev1"                    <- "lung function"
res$trait %in% "total serum ige"         <- "serum IgE"
res$trait %in% "psychosis"               <- "schizophrenia"
res$trait %in% "television viewing time" <- "sedentary behaviour"
res$trait %in% "time spent sitting"      <- "sedentary behaviour"

res <- res[order(-res$odds_diff),]


#--- save ----------------------------------------------------------------------

saveRDS(res, "ewasoverlaps_blueprint.rds")

