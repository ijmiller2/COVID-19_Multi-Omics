##### X9_KAO_Adding_Yuchens_pvalues_into_DB.R #### 

## Yuchen performed analysis on HFD for each biomolecue 
## those data are found in Rdata files in regression folder

library(DBI)
library(RSQLite)

##### read in p-values that are currently in DB ####

## Establish a connection to the DB
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 

## read Rdata files

load("E:/COVID-19_Multi-Omics/data/HFDs/df_pvalues_abundance.RData")

head(df_pvalues_abundance)
names(df_pvalues_abundance)[5] <- "formula"
df_pvalues_abundance$formula <- df_pvalues_abundance$test
df_pvalues_abundance$test<- "ANOVA"
df_pvalues_abundance$comparison <- "Hospital_free_days_45"

df_pvalues_abundance$pvalue_id <- seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_abundance), by =1 )

## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

## write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_abundance, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 

#### 2) Adding pvalues HFD ~ abundance + age less than 90

load("E:/COVID-19_Multi-Omics/data/HFDs/df_pvalues_abundance_Age_less_than_90.RData")

head(df_pvalues_abundance_Age_less_than_90)

names(df_pvalues_abundance_Age_less_than_90)[5] <- "formula"
df_pvalues_abundance_Age_less_than_90$formula <- df_pvalues_abundance_Age_less_than_90$test

df_pvalues_abundance_Age_less_than_90$test<- "ANOVA"
df_pvalues_abundance_Age_less_than_90$comparison <- "Hospital_free_days_45"

df_pvalues_abundance_Age_less_than_90$pvalue_id <- seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_abundance_Age_less_than_90), by =1 )

## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

## write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_abundance_Age_less_than_90, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 

##### 3) Adding in pvalues for HFD ~ abundance + Charlson score

load("E:/COVID-19_Multi-Omics/data/HFDs/df_pvalues_abundance_Charlson_score.RData")

head(df_pvalues_abundance_Charlson_score)

names(df_pvalues_abundance_Charlson_score)[5] <- "formula"
df_pvalues_abundance_Charlson_score$formula <- df_pvalues_abundance_Charlson_score$test
df_pvalues_abundance_Charlson_score$test<- "ANOVA"
df_pvalues_abundance_Charlson_score$comparison <- "Hospital_free_days_45"

df_pvalues_abundance_Charlson_score$pvalue_id <- seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_abundance_Charlson_score), by =1 )

## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

## write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_abundance_Charlson_score, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 

#### 4) Adding pvalues for HFD ~ abundance + ICU_1 

load("E:/COVID-19_Multi-Omics/data/HFDs/df_pvalues_abundance_ICU_1.RData")

head(df_pvalues_abundance_ICU_1)

names(df_pvalues_abundance_ICU_1)[5] <- "formula"
df_pvalues_abundance_ICU_1$formula <- df_pvalues_abundance_ICU_1$test
df_pvalues_abundance_ICU_1$test<- "ANOVA"
df_pvalues_abundance_ICU_1$comparison <- "Hospital_free_days_45"

df_pvalues_abundance_ICU_1$pvalue_id <- seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_abundance_ICU_1), by =1 )

## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

## write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_abundance_ICU_1, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 

##### 5) Adding pvalues for HFD ~ abundance + SOFA 

load("E:/COVID-19_Multi-Omics/data/HFDs/df_pvalues_abundance_SOFA.RData")

head(df_pvalues_abundance_SOFA)
names(df_pvalues_abundance_SOFA)[5] <- "formula"
df_pvalues_abundance_SOFA$formula <- df_pvalues_abundance_SOFA$test
df_pvalues_abundance_SOFA$test<- "ANOVA"
df_pvalues_abundance_SOFA$comparison <- "Hospital_free_days_45"

df_pvalues_abundance_SOFA$pvalue_id <- seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_abundance_SOFA), by =1 )

## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

## write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_abundance_SOFA, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 

#### 7) Adding pvalues for likelyhood ratio test 

load("E:/COVID-19_Multi-Omics/data/HFDs/df_pvalues_lr_Age_less_than_90.RData")
load("E:/COVID-19_Multi-Omics/data/HFDs/df_pvalues_lr_Charlson_score.RData")
load("E:/COVID-19_Multi-Omics/data/HFDs/df_pvalues_lr_ICU_1.RData")
load("E:/COVID-19_Multi-Omics/data/HFDs/df_pvalues_lr_COVID.RData")

head(df_pvalues_lr_Age_less_than_90 )

hist(df_pvalues_abundance$p_value)
hist(df_pvalues_abundance_Age_less_than_90$p_value)
hist(df_pvalues_lr_Age_less_than_90$p_value)
hist(df_pvalues_lr_Charlson_score$p_value)
hist(df_pvalues_lr_ICU_1$p_value)
hist(df_pvalues_lr_COVID$p_value)

## from Yuchen's data it looks like all these factors play a strong role in hospital free days linear regression. 
## seems like a more complex model would be appropriate 
compare_lr <- function(biomolecule_id, formula_null, formula_test, data, return = 'pvalue'){
  # comparing likelyhood ratios for models +/- fixed parameter of interest
  lm_formula_null <- lm(formula_null, data = data[data$biomolecule_id == biomolecule_id, ])
  lm_formula_test <- lm(formula_test, data = data[data$biomolecule_id == biomolecule_id, ])
  lrt <- tryCatch(anova(lm_formula_null, lm_formula_test), error = function(e) NULL)
  if (is.null(lrt)){
    lrt_lratio <- NA
    lrt_pvalue <- NA
  } else{
    lrt_lratio <- lrt$F[2]
    lrt_pvalue <- lrt$`Pr(>F)`[2]
  }
  if(return == 'pvalue') lrt_pvalue else if(return == 'both') c(lratio = lrt_lratio,pvalue = lrt_pvalue)
}


#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### Pull data from DB ####

dbListTables(con)

df_metabolites<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, metabolomics_measurements.biomolecule_id, COVID, Age_less_than_90, ICU_1, Charlson_score, SOFA, APACHEII, Hospital_free_days_45
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")

df_lipids<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, lipidomics_measurements.biomolecule_id, COVID, Age_less_than_90, ICU_1, Charlson_score, SOFA, APACHEII, Hospital_free_days_45
           FROM lipidomics_measurements
           INNER JOIN lipidomics_runs ON lipidomics_runs.replicate_id = lipidomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = lipidomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = lipidomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_proteins<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, proteomics_measurements.biomolecule_id, COVID, Age_less_than_90, ICU_1, Charlson_score, SOFA, APACHEII, Hospital_free_days_45
           FROM proteomics_measurements
           INNER JOIN proteomics_runs ON proteomics_runs.replicate_id = proteomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = proteomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = proteomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_transcripts<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, transcriptomics_measurements.biomolecule_id, COVID, Age_less_than_90, ICU_1, Charlson_score, SOFA, APACHEII, Hospital_free_days_45
           FROM transcriptomics_measurements
           INNER JOIN transcriptomics_runs ON transcriptomics_runs.replicate_id = transcriptomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = transcriptomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = transcriptomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")
dbDisconnect(con)

df <- rbind(df_metabolites, df_lipids, df_proteins, df_transcripts)

df <- df[df$sample_id != 54, ]

#### Creating dataframe to hold pvalues #######

df_pvalues <- data.frame(biomolecule_id = unique(df$biomolecule_id), test = "LR_test", comparison = "Hospital_free_days_45", formula = "Hospital_free_days_45 ~ normalized_abundance + COVID + Age_less_than_90 + Charlson_score + SOFA + APACHEII vs. Hospital_free_days_45 ~COVID + Age_less_than_90 + Charlson_score + SOFA + APACHEII ")

p_value <- apply(df_pvalues, 1, function(x)  
  compare_lr(as.numeric(x[1]), formula_null = Hospital_free_days_45 ~ COVID + Age_less_than_90 + Charlson_score + SOFA + APACHEII, 
             formula_test = Hospital_free_days_45 ~ normalized_abundance + COVID + Age_less_than_90 + Charlson_score + SOFA + APACHEII ,
             data = df, return = 'both'))

df_pvalues$effect_size <- p_value[1,]
df_pvalues$p_value <- p_value[2,]
df_pvalues$q_value <- p.adjust(df_pvalues$p_value, method = "fdr")

df_pvalues <- cbind(pvalue_id = row.names(df_pvalues), df_pvalues)

df_pvalues$pvalue_id <- seq(nrow(pvalues)+1, length.out = nrow(df_pvalues), by =1 )

hist(df_pvalues$p_value)

## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

## write table to DB 

dbWriteTable(con, "pvalues", df_pvalues, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 
