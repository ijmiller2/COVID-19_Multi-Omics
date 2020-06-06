####### 11_KAO_Looking_at_effect_of_DM_status.R ####### 

library(DBI)
library(RSQLite)

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

lm_coefficients <- function(biomolecule_id, formula_test, data, return_first = T, return_category = NULL){
  lm_formula_test <- lm(formula_test, data = data[data$biomolecule_id == biomolecule_id, ])
  if(return_first == T){
    c(coefficent =  summary(lm_formula_test)$coefficients[2,3], pvalue = summary(lm_formula_test)$coefficients[2,4])
  } else if (!is.null(return_category) & return_category %in% row.names(summary(lm_formula_test)$coefficients)){
    index <- match(return_category, row.names(summary(lm_formula_test)$coefficients))
    c(coeficient =  summary(lm_formula_test)$coefficients[index,3], pvalue = summary(lm_formula_test)$coefficients[index,4])
  }
}


anova_effect <- function(biomolecule_id, formula_test, data, return_first = T, return_index = NULL){
  lm_formula_test <- lm(formula_test, data = data[data$biomolecule_id == biomolecule_id, ])
  if(return_first == T){
    c(effect_size =  anova(lm_formula_test)$F[1], pvalue = anova(lm_formula_test)$P[1])==
  } else if (!is.null(return_index)){
    c(effect_size =  anova(lm_formula_test)$F[index], pvalue = anova(lm_formula_test)$P[index])
  }
}

###### Applying this function across all features ########

#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### Pull data from DB ####

dbListTables(con)

df_metabolites<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, metabolomics_measurements.biomolecule_id, COVID, Hospital_free_days_45, Age_less_than_90, Gender, ICU_1, DM, Charlson_score, SOFA
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")

df_lipids<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, lipidomics_measurements.biomolecule_id, COVID, Hospital_free_days_45, Age_less_than_90, Gender, ICU_1, DM, Charlson_score, SOFA
           FROM lipidomics_measurements
           INNER JOIN lipidomics_runs ON lipidomics_runs.replicate_id = lipidomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = lipidomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = lipidomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_proteins<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, proteomics_measurements.biomolecule_id, COVID, Hospital_free_days_45, Age_less_than_90, Gender, ICU_1, DM, Charlson_score, SOFA
           FROM proteomics_measurements
           INNER JOIN proteomics_runs ON proteomics_runs.replicate_id = proteomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = proteomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = proteomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")
dbDisconnect(con)


df <- rbind(df_metabolites, df_lipids, df_proteins)

df <- df[df$sample_id != 54, ]

df$SOFA <- as.numeric(df$SOFA)

###### HOstpital free days and other covariates ##### 

# simple model 
fit1 <- lm(Hospital_free_days_45 ~ COVID, df[df$biomolecule_id == 1, ])
anova(fit1)

fit2 <- lm(Hospital_free_days_45 ~ COVID * DM, df[df$biomolecule_id == 1, ])
anova(fit2)

anova(fit1, fit2)

# more complex model 
fit3 <- lm(Hospital_free_days_45 ~ COVID + ICU_1 + Charlson_score + Age_less_than_90 , data = df[df$biomolecule_id == 1, ])
anova(fit3)

fit4 <- lm(Hospital_free_days_45 ~ COVID + Charlson_score + Age_less_than_90 , data = df[df$biomolecule_id == 1, ])
anova(fit4)

###### fitting Simple model with biomolecules ###### 
df_pvalues <- df_pvalues <- data.frame(biomolecule_id = unique(df$biomolecule_id), test = "LR_test", comparison = "Biomolecule_id_effect_on_HFD", confounders = "COVID")

pvalues_simple <- apply(df_pvalues, 1, function(x)  
  compare_lr(as.numeric(x[1]), formula_null = Hospital_free_days_45 ~ normalized_abundance + COVID, 
             formula_test = Hospital_free_days_45 ~ COVID,
             data = df, return = 'pvalue'))

hist(pvalues_simple, breaks = 100, main = "Hospital free days 45 ~ biomolecule + COVID")


pvalues_simple_DM <- apply(df_pvalues, 1, function(x)  
  compare_lr(as.numeric(x[1]), formula_null = Hospital_free_days_45 ~ normalized_abundance + COVID + DM, 
             formula_test = Hospital_free_days_45 ~ COVID + DM,
             data = df, return = 'pvalue'))

hist(pvalues_simple_DM, breaks = 100, main = "Hospital free days 45 ~ biomolecule + COVID + DM")

df_pvalues$biomolecule_id[(pvalues_simple < 0.05) != (pvalues_simple_DM < 0.05 )] ## first molecule looks like glucose (Sugar rt 16.65 on GC)
pvalues_simple[(pvalues_simple < 0.05) != (pvalues_simple_DM < 0.05 )] ## all these p values are borderline ~ 0.05 
pvalues_simple_DM[(pvalues_simple < 0.05) != (pvalues_simple_DM < 0.05 )] ## accounting for DM only changes them slightly (still borderline 0.05)

pvalues_simple_DM_effect <- apply(df_pvalues, 1, function(x)  
  compare_lr(as.numeric(x[1]), formula_null = Hospital_free_days_45 ~ normalized_abundance + COVID , 
             formula_test = Hospital_free_days_45 ~ normalized_abundance + COVID + DM,
             data = df, return = 'pvalue'))


hist(pvalues_simple_DM_effect, breaks = 100, main = "Hospital free days 45 ~ biomolecule + COVID + DM \nvs.Hospital free days 45 ~ biomolecule + COVID  ")

