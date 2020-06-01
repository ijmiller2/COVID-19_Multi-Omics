###### X7_KAO_updating_metadata_biomolecule_id.R ##### 

## Updating metadata for lipids because biomolecule_id in metadata is incorrect. 


library(DBI)
library(RSQLite)

#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### Pull data from DB ####

##### 
df<- dbGetQuery(con, "SELECT unique_identifier, normalized_abundance, lipidomics_measurements.biomolecule_id, batch 
           FROM lipidomics_measurements
           INNER JOIN lipidomics_runs ON lipidomics_runs.replicate_id = lipidomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = lipidomics_runs.rawfile_id
           WHERE keep = 1
           ")

biomolecule_df <- dbGetQuery(con, "SELECT *
                             FROM biomolecules
                             INNER JOIN metadata on metadata.biomolecule_id = biomolecules.biomolecule_id
                             WHERE omics_id = 2", stringsAsFactors = F)

metadata_df <- dbGetQuery(con, "SELECT *
                       FROM metadata
                       WHERE metadata_type = 'Lipid Class'")
dbDisconnect(con)

##### Load in data w/ slope of control dilution and whether the feature was found in control ##### 

df_slope <- read.csv("P:/All_20200428_COVID_plasma_multiomics/Lipidomics/Lipidomics_quant_results/Final_Results_Slopes_KAO.csv", stringsAsFactors = F)

##### Converting to wide format ######

# being explicit about items in the data frame, in case db changes 
biomolecule_df_2 <- data.frame(metdata_id = as.numeric(biomolecule_df$metadata_id),
                               biomolecule_id = rep(unique(df$biomolecule_id), each = length(unique(biomolecule_df$metadata_type) )), 
                               standardized_name = as.character(biomolecule_df$standardized_name),
                               metadata_type = as.character(biomolecule_df$metadata_type), 
                               metadata_value = as.character(biomolecule_df$metadata_value),
                               stringsAsFactors = F)

bio_df_wide <- reshape(biomolecule_df_2[,-1], v.names = "metadata_value", timevar = "metadata_type", idvar = c("biomolecule_id"), direction = "wide")
row.names(bio_df_wide) <- bio_df_wide$biomolecule_id

# checking that these dimentions match with the df_slope dim 
dim(bio_df_wide) ## incorrect;  7235 8
dim(df_slope) # 7235 11
length(unique(df$biomolecule_id)) #7235
dim(metadata_df) # 7235 (duplicates here too, should)

##### df_slope table contains correct metadata, need to update metadata table #### 

# add biomolecule id to df_slope
df_slope$biomolecule_id  <- df$biomolecule_id[1:length(unique(df$biomolecule_id))]

# checking these tables match 
all(df_slope$biomolecule_id == bio_df_wide$biomolecule_id)

# checking that RT matches; not true, but only b/c of rounding in df_slope
table(df_slope$Retention_Time__min_ == bio_df_wide$`metadata_value.Retention Time (min)`)
bio_df_wide[df_slope$Retention_Time__min_ != bio_df_wide$`metadata_value.Retention Time (min)`,]
df_slope[df_slope$Retention_Time__min_ != bio_df_wide$`metadata_value.Retention Time (min)`,]

###### Update metadata table in db #######

df_replace <- data.frame(metadata_id = biomolecule_df_2$metdata_id, 
                         biomolecule_id = biomolecule_df_2$biomolecule_id)

df_replace <- df_replace[biomolecule_df$biomolecule_id != biomolecule_df_2$biomolecule_id, ]

##### R Function to replace values ##### 
replaceValues <- function(x) dbExecute(con, paste("UPDATE metadata
          SET biomolecule_id = ",unlist(x[2])," WHERE metadata_id = ", unlist(x[1]), sep = ""))

##### Establish connection to DB #######

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

##### Iterate over data frame to update normalized_abundance values in sqlite db ######

apply(df_replace, 1, replaceValues)

##### Check that values match ##### 
metadata_df <- dbGetQuery(con, "SELECT *
                       FROM metadata
                       WHERE metadata_type = 'Lipid Class'")

dbDisconnect(con)

length(metadata_df$biomolecule_id) #7235
length(unique(metadata_df$biomolecule_id)) #7235

#####

