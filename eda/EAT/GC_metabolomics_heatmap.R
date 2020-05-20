#### 01_EAT_exploring_KAO_04_GCdata.R ##### 

## Updates to the database help with feature filtering of the 
## GC metabolites. See "X4_KAO_Updating_biomolecules_keep_column
## _for_GC_metabolites.R" and "03_KAO_Exploring_GC_feature_quality.R"
## 


## Access DB SQlite where GC data lives to use in pheatmap. 
## The code was taken from KAO 

library(DBI)
library(RSQLite)
library(pheatmap)
library(plyr)

##### Establish a connection with DB #####

con <- dbConnect(SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

##### Pull data from DB #####

dbListTables(con)
# [1] "biomolecules"                  "deidentified_patient_metadata" "lipidomics_measurements"       "lipidomics_runs"              
# [5] "metabolomics_measurements"     "metabolomics_runs"             "metadata"                      "omes"                         
# [9] "patient_metadata"              "patient_samples"               "proteomics_measurements"       "proteomics_runs"              
# [13] "rawfiles"                      "sqlite_sequence"              

dbReadTable(con, "biomolecules")[0,]
dbReadTable(con, "metadata")[1:10,]
dbReadTable(con, "metabolomics_measurements")[0,]
dbReadTable(con, "metabolomics_runs")[0,]
dbReadTable(con, "biomolecules")[0,]

##### 
df<- dbGetQuery(con, "SELECT unique_identifier, normalized_abundance, metabolomics_measurements.biomolecule_id, batch 
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1 
           AND ome_id = 3 
           AND biomolecules.keep = '1'
           ")
dbDisconnect(con)


df_wide <- reshape(df, timevar = "biomolecule_id", v.names = "normalized_abundance",
                   idvar = "unique_identifier", direction = "wide" )

meta <- as.data.frame(stringr::str_split_fixed(df_wide$unique_identifier,"_",5))
meta$V5 <- meta$V4
meta$V6 <- df_wide$batch
colnames(meta) <- c("Date", "Sample.Type", "Disease.State", "Patient.ID.unique","Patient.ID.group","Batch")

meta$Patient.ID.unique <- make.unique(as.character(meta$Patient.ID.unique))
meta$Disease.State[which(meta$Disease.State == "Conrol")] <- "Control"
meta$Disease.State <- factor(meta$Disease.State,levels = c("Control","COVID","NONCOVID"))


metabolomics <- df_wide[,c(3:ncol(df_wide))]
rownames(metabolomics) <- meta$Patient.ID.unique
colnames(metabolomics) <- sub(".*abundance.","",colnames(metabolomics))
pheatmap::pheatmap(metabolomics)
