#### 13_KAO_calculating_FC_for_COVID.R ##### 

## data table with effect size for COVID #### 

library(DBI)
library(RSQLite)


lm_coefficients <- function(biomolecule_id, formula_test, data, return_first = T, return_category = NULL){
  lm_formula_test <- lm(formula_test, data = data[data$biomolecule_id == biomolecule_id, ])
  if(return_first == T){
    c(coeficient =  summary(lm_formula_test)$coefficients[2,3], pvalue = summary(lm_formula_test)$coefficients[2,4])
  } else if (!is.null(return_category) & return_category %in% row.names(summary(lm_formula_test)$coefficients)){
    index <- match(return_category, row.names(summary(lm_formula_test)$coefficients))
    c(coeficient =  summary(lm_formula_test)$coefficients[index,3], pvalue = summary(lm_formula_test)$coefficients[index,4])
  }
}

###### Applying this function across all features ########

#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### Pull data from DB ####

dbListTables(con)

df_metabolites<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, metabolomics_measurements.biomolecule_id, COVID, Age_less_than_90, Gender, ICU_1, DM, Charlson_score, SOFA
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")

df_lipids<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, lipidomics_measurements.biomolecule_id, COVID, Age_less_than_90, Gender, ICU_1, DM, Charlson_score, SOFA
           FROM lipidomics_measurements
           INNER JOIN lipidomics_runs ON lipidomics_runs.replicate_id = lipidomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = lipidomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = lipidomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_proteins<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, proteomics_measurements.biomolecule_id, COVID, Age_less_than_90, Gender, ICU_1, DM, Charlson_score, SOFA
           FROM proteomics_measurements
           INNER JOIN proteomics_runs ON proteomics_runs.replicate_id = proteomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = proteomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = proteomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_transcripts<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, transcriptomics_measurements.biomolecule_id, COVID, Age_less_than_90, Gender, ICU_1, DM, Charlson_score, SOFA
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

df_COVID <- data.frame(biomolecule_id = unique(df$biomolecule_id), test = "lm_coefficients", comparison = "COVID_vs_NONCOVID", formula = "normalized_abundance ~ COVID")

COVID_effects <- apply(df_COVID, 1, function(x)  
  lm_coefficients(as.numeric(x[1]), 
             formula_test = normalized_abundance ~ COVID,
             data = df))

df_COVID$effect_size <- COVID_effects[1,]

plot(df_COVID$effect_size, -log(pvalues_1$p_value))

means<- aggregate(df$normalized_abundance, by = list(df$COVID, df$biomolecule_id), mean)
head(means)
FC <- data.frame(FC = means$x[as.integer(row.names(means)) %% 2 == 0] - means$x[as.integer(row.names(means)) %% 2 == 1] )
row.names(FC) <- means$Group.2[as.integer(row.names(means)) %% 2 == 0]

m <- merge(df_COVID, FC, by.x = 'biomolecule_id', by.y = 0)
plot(m$effect_size, m$FC)

write.csv(m, file = "data/COVID_fc_by_biomolecule_ID.csv")

