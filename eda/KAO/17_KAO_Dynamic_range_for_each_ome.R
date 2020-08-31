##### Dynamic range for each ome ##### 

## These data pull infromation from pvalues table and provides numbers for the text about figure 2. 

library(DBI)
library(RSQLite)

source("eda/KAO/0_pathway_toolkit.R")
library(pheatmap)
library(ellipse)
library(scales)

colors <- read.csv("reference/color_palette.txt",  stringsAsFactors = F)[,2]
palette(colors)

### Fig 2A - PCA #### 

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")


df_metabolites<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, metabolomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")

df_lipids<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, lipidomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45
           FROM lipidomics_measurements
           INNER JOIN lipidomics_runs ON lipidomics_runs.replicate_id = lipidomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = lipidomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = lipidomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_proteins<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, proteomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45
           FROM proteomics_measurements
           INNER JOIN proteomics_runs ON proteomics_runs.replicate_id = proteomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = proteomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = proteomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_transcripts<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, transcriptomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45
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
names(df)

#### Creating a wide-format data frame to facilitate PCA #####

df_wide <- reshape(df, timevar = "biomolecule_id", v.names = "normalized_abundance",
                   idvar = "sample_id", direction = "wide" )

names(df_wide)
df_wide_all <- df_wide[!(df_wide$sample_id == 88 | df_wide$sample_id == 81 | df_wide$sample_id == 129), ]

df_wide_exprs <- df_wide_all[,-c(1:4)]


boxplot(colMeans(df_wide_exprs[, names(df_wide_exprs) %in% m2$biomolecule_id[m2$omics_id.x == 1]]),
        colMeans(df_wide_exprs[, names(df_wide_exprs) %in% m2$biomolecule_id[m2$omics_id.x == 2]]),
colMeans(df_wide_exprs[, names(df_wide_exprs) %in% m2$biomolecule_id[m2$omics_id.x == 3]]),
colMeans(df_wide_exprs[, names(df_wide_exprs) %in% m2$biomolecule_id[m2$omics_id.x == 5]]),
names = c("Proteins", "Lipids", " Metabolites", "Transcipts"))

boxplot(df_wide_all)
