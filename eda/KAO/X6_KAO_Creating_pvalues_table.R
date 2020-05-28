##### X6_KAO_creating_pvalue_table.R ######

## The goal of this script is to develop a function to produce model comparisons 
## between models with COVID vs. not and output likelyhood ratio pvalues. This 
## approach allows for inclusion of confounders.
## 
## These data will be stored in a table in the database called pvalues
## 

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
  if(return == 'pvalue') lrt_pvalue else lrt_lratio
}

###### Applying this function across all features ########

#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### Pull data from DB ####

dbListTables(con)


##### 
df_metabolites<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, metabolomics_measurements.biomolecule_id, COVID, Age_less_than_90, Gender, ICU_1
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")

df_lipids<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, lipidomics_measurements.biomolecule_id, COVID, Age_less_than_90, Gender, ICU_1
           FROM lipidomics_measurements
           INNER JOIN lipidomics_runs ON lipidomics_runs.replicate_id = lipidomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = lipidomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = lipidomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_proteins<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, proteomics_measurements.biomolecule_id, COVID, Age_less_than_90, Gender, ICU_1
           FROM proteomics_measurements
           INNER JOIN proteomics_runs ON proteomics_runs.replicate_id = proteomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = proteomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = proteomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")
dbDisconnect(con)

#### Creating dataframe to hold pvalues #######

df <- rbind(df_metabolites,df_lipids, df_proteins)

df_pvalues <- data.frame(biomolecule_id = unique(df$biomolecule_id), test = "LR_test", comparison = "COVID_vs_NONCOVID", confounders = "ICU_1;Gender;Age_less_than_90")

df_pvalues$p_value <- apply(df_pvalues, 1, function(x)  
            compare_lr(as.numeric(x[1]), formula_null = normalized_abundance ~ ICU_1 + Gender + Age_less_than_90, 
             formula_test = normalized_abundance ~ COVID * ICU_1 + Gender + Age_less_than_90,
             data = df, return = 'pvalue'))

df_pvalues$q_value <- p.adjust(df_pvalues$p_value, method = "fdr")

df_pvalues <- cbind(pvalue_id = row.names(df_pvalues), df_pvalues)

#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### write table to DB ####

dbWriteTable(con, "pvalues", df_pvalues)

# check
dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 

