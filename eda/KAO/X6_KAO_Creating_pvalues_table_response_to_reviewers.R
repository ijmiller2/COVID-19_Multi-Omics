##### X6_KAO_creating_pvalue_table.R ######

## The goal of this script is to develop a function to produce model comparisons 
## between models with COVID vs. not and output likelyhood ratio pvalues. This 
## approach allows for inclusion of confounders.
## 
## These data will be stored in a table in the database called pvalues
## 
## 20200601 - KAO edit due to earlier error in sample IDs - this new calculation
## overwrites the old table. 
## 20200828 - KAO edit to recalculate with robust linear regression model. and also calculate severity based on WHO.

library(DBI)
library(RSQLite)
library(MASS)
library(lmtest)

##### Function to generate a likelyhood ratio based p-value #####

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

compare_lr_robust <- function(biomolecule_id, formula_null, formula_test, data, return = 'pvalue'){
  # comparing likelyhood ratios for models +/- fixed parameter of interest
  set.seed(1234)
  lm_formula_null <- rlm(formula_null, data = data[data$biomolecule_id == biomolecule_id, ])
  set.seed(1234)
  lm_formula_test <- rlm(formula_test, data = data[data$biomolecule_id == biomolecule_id, ])
  lrt <- tryCatch(lrtest(lm_formula_null, lm_formula_test), error = function(e) NULL)
  if (is.null(lrt)){
    lrt_lratio <- NA
    lrt_pvalue <- NA
  } else{
    lrt_lratio <- lrt$Chisq[2]
    lrt_pvalue <- lrt$`Pr(>Chisq)`[2]
  }
  if(return == 'pvalue') lrt_pvalue else if(return == 'both') c(lratio = lrt_lratio,pvalue = lrt_pvalue)
}

# this function is currently usign SS type III calculation. But probably we should be doing SS type I. Will edit in future. 
# lm_coefficients <- function(biomolecule_id, formula_test, data, return_first = T, return_category = NULL){
#   lm_formula_test <- lm(formula_test, data = data[data$biomolecule_id == biomolecule_id, ])
#   if(return_first == T){
#     c(coeficient =  summary(lm_formula_test)$coefficients[2,3], pvalue = summary(lm_formula_test)$coefficients[2,4])
#   } else if (!is.null(return_category) & return_category %in% row.names(summary(lm_formula_test)$coefficients)){
#     index <- match(return_category, row.names(summary(lm_formula_test)$coefficients))
#     c(coeficient =  summary(lm_formula_test)$coefficients[index,3], pvalue = summary(lm_formula_test)$coefficients[index,4])
#   }
# }

###### Applying this function across all features ########

#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### Pull data from DB ####

dbListTables(con)

df_metabolites<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, metabolomics_measurements.biomolecule_id, COVID, Age_less_than_90, Gender, ICU_1, Charlson_score, SOFA, WHO_ordinal_at_day_28, Hospital_free_days_45
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")

df_lipids<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, lipidomics_measurements.biomolecule_id, COVID, Age_less_than_90, Gender, ICU_1, Charlson_score, SOFA, WHO_ordinal_at_day_28, Hospital_free_days_45
           FROM lipidomics_measurements
           INNER JOIN lipidomics_runs ON lipidomics_runs.replicate_id = lipidomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = lipidomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = lipidomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_proteins<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, proteomics_measurements.biomolecule_id, COVID, Age_less_than_90, Gender, ICU_1, Charlson_score, SOFA, WHO_ordinal_at_day_28, Hospital_free_days_45
           FROM proteomics_measurements
           INNER JOIN proteomics_runs ON proteomics_runs.replicate_id = proteomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = proteomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = proteomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_transcripts<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, transcriptomics_measurements.biomolecule_id, COVID, Age_less_than_90, Gender, ICU_1, Charlson_score, SOFA, WHO_ordinal_at_day_28, Hospital_free_days_45
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

# df_pvalues <- data.frame(biomolecule_id = unique(df$biomolecule_id), test = "LR_test", comparison = "COVID_vs_NONCOVID", formula = "normalized_abundance ~ COVID * ICU_1 + Gender + Age_less_than_90 vs. normalized_abundance ~ ICU_1 + Gender + Age_less_than_90")
# 
# p_value <- apply(df_pvalues, 1, function(x)  
#             compare_lr(as.numeric(x[1]), formula_null = normalized_abundance ~ ICU_1 + Gender + Age_less_than_90, 
#              formula_test = normalized_abundance ~ COVID * ICU_1 + Gender + Age_less_than_90,
#              data = df, return = 'both'))
# 
# df_pvalues$effect_size <- p_value[1,]
# df_pvalues$p_value <- p_value[2,]
# df_pvalues$q_value <- p.adjust(df_pvalues$p_value, method = "fdr")
# 
# df_pvalues <- cbind(pvalue_id = row.names(df_pvalues), df_pvalues)
# 
# hist(df_pvalues$p_value)
# 
# #### Establish a connection to the DB
# con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")
# 
# #### write table to DB 
# 
# dbWriteTable(con, "pvalues", df_pvalues, overwrite = T)
# 
# # check
# pvalues <- dbReadTable(con, "pvalues")
# 
# # disconnect
# dbDisconnect(con) 
# 
# ##### P-values for gender ###### 
# df_gender <- df[df$Gender != "", ]
# 
# df_pvalues_gender <- data.frame(biomolecule_id = unique(df$biomolecule_id), test = "LR_test", comparison = "GENDER", formula = "normalized_abundance ~ COVID + ICU_1 + Gender + Age_less_than_90 vs. normalized_abundance ~ COVID + ICU_1 + Age_less_than_90")
# 
# p_value <- apply(df_pvalues_gender, 1, function(x)  
# compare_lr(as.numeric(x[1]), formula_null = normalized_abundance ~ COVID + ICU_1 + Age_less_than_90, 
#              formula_test = normalized_abundance ~ COVID + ICU_1 + Gender + Age_less_than_90,
#              data = df_gender, return = 'both'))
# df_pvalues_gender$effect_size <- p_value[1,]
# df_pvalues_gender$p_value <- p_value[2,]
# 
# hist(df_pvalues_gender$p_value, breaks = 100, main = 'Histogram of pvalues +/- gender')
# 
# df_pvalues_gender$q_value <- p.adjust(df_pvalues_gender$p_value, method = "fdr")
# 
# df_pvalues_gender <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_gender), by =1 ), df_pvalues_gender)
# 
# ## Append gender Pvalues to db 
# ## Establish a connection to the DB 
# con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")
# 
# #### write table to DB 
# 
# dbWriteTable(con, "pvalues", df_pvalues_gender, append = T)
# 
# # check
# pvalues <- dbReadTable(con, "pvalues")
# 
# # disconnect
# dbDisconnect(con) 
# 
# ####### P-values w/ age ######
# 
# df_pvalues_age <- data.frame(biomolecule_id = unique(df$biomolecule_id), test = "LR_test", comparison = "Age_less_than_90", formula = "normalized_abundance ~ COVID + ICU_1 + Gender + Age_less_than_90 vs. normalized_abundance ~ COVID + ICU_1 + Gender")
# 
# p_value <- apply(df_pvalues_age, 1, function(x)  
#   compare_lr(as.numeric(x[1]), formula_null = normalized_abundance ~ COVID + ICU_1 + Gender, 
#              formula_test = normalized_abundance ~ COVID + ICU_1 + Gender + Age_less_than_90,
#              data = df, return = 'both'))
# 
# df_pvalues_age$effect_size <- p_value[1,]
# df_pvalues_age$p_value <- p_value[2,]
# 
# hist(df_pvalues_age$p_value, breaks = 100, main = 'Histogram of pvalues +/- age')
# 
# 
# df_pvalues_age$q_value <- p.adjust(df_pvalues_age$p_value, method = "fdr")
# 
# df_pvalues_age <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_age), by =1 ), df_pvalues_age)
# 
# ## Append age Pvalues to db 
# ## Establish a connection to the DB 
# con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")
# 
# #### write table to DB 
# 
# dbWriteTable(con, "pvalues", df_pvalues_age, append = T)
# 
# # check
# pvalues <- dbReadTable(con, "pvalues")
# 
# # disconnect
# dbDisconnect(con) 
# 
# ####### P-values w/ ICU status ######
# 
# df_pvalues_icu <- data.frame(biomolecule_id = unique(df$biomolecule_id), test = "LR_test", comparison = "ICU", formula = "normalized_abundance ~ COVID + ICU_1 + Gender + Age_less_than_90 vs. normalized_abundance ~ COVID + Gender + Age_less_than_90")
# 
# p_value <- apply(df_pvalues_icu, 1, function(x)  
#   compare_lr(as.numeric(x[1]), formula_null = normalized_abundance ~ COVID + Gender + Age_less_than_90, 
#              formula_test = normalized_abundance ~ COVID + ICU_1 + Gender + Age_less_than_90,
#              data = df, return = 'both'))
# df_pvalues_icu$effect_size <- p_value[1,]
# df_pvalues_icu$p_value <- p_value[2,]
# 
# hist(df_pvalues_icu$p_value, breaks = 100, main = 'Histogram of pvalues +/- ICU')
# 
# 
# df_pvalues_icu$q_value <- p.adjust(df_pvalues_icu$p_value, method = "fdr")
# 
# df_pvalues_icu <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_icu), by =1 ), df_pvalues_icu)
# 
# ## Append icu Pvalues to db 
# ## Establish a connection to the DB 
# con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")
# 
# #### write table to DB 
# 
# dbWriteTable(con, "pvalues", df_pvalues_icu, append = T)
# 
# # check
# pvalues <- dbReadTable(con, "pvalues")
# 
# # disconnect
# dbDisconnect(con) 
# 
# 
# ####### Pvalues w/ COVID ICU interaction #####
# 
# df_pvalues_interaction <- data.frame(biomolecule_id = unique(df$biomolecule_id), test = "LR_test", comparison = "COVID ICU interaction", formula = "normalized_abundance ~ COVID * ICU_1 + Gender + Age_less_than_90 vs. normalized_abundance ~ COVID + ICU_1 + Gender + Age_less_than_90")
# 
# p_value <- apply(df_pvalues_interaction, 1, function(x)  
#   compare_lr(as.numeric(x[1]), formula_null = normalized_abundance ~ COVID + ICU_1 + Gender + Age_less_than_90, 
#              formula_test = normalized_abundance ~ COVID * ICU_1 + Gender + Age_less_than_90,
#              data = df, return = 'both'))
# df_pvalues_interaction$effect_size <- p_value[1,]
# df_pvalues_interaction$p_value <- p_value[2,]
# 
# hist(df_pvalues_interaction$p_value, breaks = 100, main = 'Histogram of pvalues +/- COVID ICU interaction')
# 
# df_pvalues_interaction$q_value <- p.adjust(df_pvalues_interaction$p_value, method = "fdr")
# 
# df_pvalues_interaction <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_interaction), by =1 ), df_pvalues_interaction)
# 
# ## Append interaction Pvalues to db 
# ## Establish a connection to the DB 
# con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")
# 
# #### write table to DB 
# 
# dbWriteTable(con, "pvalues", df_pvalues_interaction, append = T)
# 
# # check
# pvalues <- dbReadTable(con, "pvalues")
# 
# # disconnect
# dbDisconnect(con) 

####### P-values w/ COVID status using robust lm (rlm) ######
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con)


df_pvalues_covid_robust <- data.frame(biomolecule_id = unique(df$biomolecule_id), test = "LR_test_robust", comparison = "COVID", formula = "normalized_abundance ~ COVID + ICU_1 + Gender + Age_less_than_90 vs. normalized_abundance ~ ICU_1 + Gender + Age_less_than_90")

p_value <- apply(df_pvalues_covid_robust, 1, function(x)  
  compare_lr_robust(as.numeric(x[1]), formula_null = normalized_abundance ~ ICU_1 + Gender + Age_less_than_90, 
             formula_test = normalized_abundance ~ COVID + ICU_1 + Gender + Age_less_than_90,
             data = df, return = 'both'))
df_pvalues_covid_robust$effect_size <- p_value[1,]
df_pvalues_covid_robust$p_value <- p_value[2,]

hist(df_pvalues_covid_robust$p_value, breaks = 100, main = 'Histogram of pvalues +/- COVID')


df_pvalues_covid_robust$q_value <- p.adjust(df_pvalues_covid_robust$p_value, method = "fdr")

df_pvalues_covid_robust <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_covid_robust), by =1 ), df_pvalues_covid_robust)

## Append COVID Pvalues to db 
## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_covid_robust, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 

#### Adding HFD pvalues using robust lm #####

df_pvalues_HFD <- data.frame(biomolecule_id = unique(df$biomolecule_id), test = "LR_test_robust", comparison = "Hospital_free_days_45", formula = "Hospital_free_days_45 ~ normalized_abundance + Age_less_than_90 + Gender vs. Hospital_free_days_45 ~ Age_less_than_90 + Gender")

p_value <- apply(df_pvalues_HFD, 1, function(x)
  compare_lr_robust(as.numeric(x[1]), formula_null = Hospital_free_days_45 ~ Age_less_than_90 + Gender,
             formula_test = Hospital_free_days_45 ~ normalized_abundance + Age_less_than_90 + Gender,
             data = df, return = 'both'))


df_pvalues_HFD$effect_size <- p_value[1,]
df_pvalues_HFD$p_value <- p_value[2,]
hist(df_pvalues_HFD$p_value)

df_pvalues_HFD$q_value <- p.adjust(df_pvalues_HFD$p_value, method = "fdr")

df_pvalues_HFD <- cbind(pvalue_id = seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_HFD), by =1 ), df_pvalues_HFD)


#
# ## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

## write table to DB

dbWriteTable(con, "pvalues", df_pvalues_HFD, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con)

