##### X4_KAO_updating_biomolecules_keep_column_GC_metabolites.R" ######

## Some GC-metabolites should be excluded because they are 
## duplicates and the overall RSDs are over 0.3. See script
## eda/KAO/03_KAO_Exploring_GC_feature_quality.R
## 
## Biomolcules 'keep' column is default to 1; for biomoleues 
## which should not be kept, will update column to "0; reason for
## exclusion 1; reason for exclusion 2; etc." 
##
## CAUTION: this script modifies db and iterates through script
## do not run this script again without checking biomolecules keep. 

library(DBI)
library(RSQLite)

#### Establish a connection to the DB #####

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

##### checking table structure #### 

dbListTables(con)

##### Extracting the biomolecule identifier and relavant info from sqlite db ##### 

df<- dbGetQuery(con, "SELECT unique_identifier, normalized_abundance, metabolomics_measurements.biomolecule_id, batch 
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN metadata on metadata.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE keep = 1 
           AND ome_id = 3 
           ")


biomolecule_df <- dbGetQuery(con, "SELECT *
                             FROM biomolecules
                             INNER JOIN metadata on metadata.biomolecule_id = biomolecules.biomolecule_id
                             WHERE omics_id = 3
                             AND (metadata_type = 'Remove_duplicate') ", stringsAsFactors = F)

dbDisconnect(con)

##### For biomolecules which are duplicated, update keep for duplicates to "0; is_duplicate" #######

# making a smaller data frame which only contains biomoleucles that are duplicates
duplicate <- biomolecule_df$metadata_value == "T"
biomolecule_subset <- data.frame(biomolecule_id = biomolecule_df$biomolecule_id[duplicate], keep = biomolecule_df$keep[duplicate], stringsAsFactors = F)

# modify the keep string; if 1 turn to 0, then append "is_duplicate"
biomolecule_subset$keep[biomolecule_subset$keep == "1"] <- sub("1", "0", biomolecule_subset$keep[biomolecule_subset$keep == "1"])
biomolecule_subset$keep <- paste(biomolecule_subset$keep, "is_duplicate", sep = ";")        

# funtion to replace biomolecules keep values
replaceValues <- function(x)  dbExecute( con, paste("UPDATE biomolecules
          SET keep = '",unlist(x[2]),"' WHERE biomolecule_id = ", x[1], sep = "")) 

##### Iterate over smaller data frame to update biomolecule keep values in sqlite db ######

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#dbExecute(con, "UPDATE biomolecules SET keep = '1'")                        

apply(biomolecule_subset, 1, replaceValues)

## confirm 
biomolecule_df <- dbGetQuery(con, "SELECT *
          FROM biomolecules
          ")

dbDisconnect(con)

##### determine if intrabatch QC csv are > 0.3 and update keep to 0 #####

## FUNCTION to calculate RSD 
rsd <- function(x) sd(2^x)/mean(2^x)

# select for controls
controls <- grepl("Control", df$unique_identifier)

# calculate across batch QC rsds 
rsd_controls <- aggregate(df$normalized_abundance[controls], 
                          by = list(biomolecule_id = df$biomolecule_id[controls]),
                          rsd)

# create a subset of molecules with > 0.3 RSD 
rsd_controls_subset <- rsd_controls[rsd_controls$x > 0.3,]
names(rsd_controls_subset) <- c("biomolecule_id", "keep")
rsd_controls_subset$keep <- biomolecule_df$keep[match(rsd_controls_subset$biomolecule_id, biomolecule_df$biomolecule_id)]

# modify the keep string; if 1 turn to 0, then append "QC_interbatch_RSD_over30perc"
rsd_controls_subset$keep[rsd_controls_subset$keep == "1"] <- sub("1", "0", rsd_controls_subset$keep[rsd_controls_subset$keep == "1"])
rsd_controls_subset$keep <- paste(rsd_controls_subset$keep, "QC_interbatch_RSD_over30perc", sep = ";")        

##### Iterate over smaller data frame to update biomolecule keep values in sqlite db ######

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#dbExecute(con, "UPDATE biomolecules SET keep = '1'")                        

apply(rsd_controls_subset, 1, replaceValues)

## confirm 
biomolecule_df <- dbGetQuery(con, "SELECT *
          FROM biomolecules
          ")

dbDisconnect(con)

##### update keep to 0 if 2 or more intra batch QC RSDs are greater than 0.3 ######


# calculate rsd for controls by biomolecule and batch
rsd_controls_by_batch <- aggregate(df$normalized_abundance[controls], 
                                   by = list(batch = df$batch[controls], biomolecule_id = df$biomolecule_id[controls]), 
                                   rsd)

# reshape by batch into wide format 
rsd_controls_by_batch_wide <- reshape(rsd_controls_by_batch, v.names = "x", timevar = "batch", idvar = c("biomolecule_id"), direction = "wide")
names(rsd_controls_by_batch_wide) <- sub("x.", "rsd_batch_", names(rsd_controls_by_batch_wide))

# create subset of features which need keep updated 
by_batch_subset <- data.frame( biomolecule_id = rsd_controls_by_batch_wide$biomolecule_id[rowSums(rsd_controls_by_batch_wide[,-1]>0.3)>1],keep = biomolecule_df$keep[match(rsd_controls_by_batch_wide$biomolecule_id[rowSums(rsd_controls_by_batch_wide[,-1]>0.3)>1], biomolecule_df$biomolecule_id)], stringsAsFactors = F)

# modify the keep string; if 1 turn to 0, then append "QC_2_intrabatch_RSD_over30perc"
by_batch_subset$keep[by_batch_subset$keep == "1"] <- sub("1", "0", by_batch_subset$keep[by_batch_subset$keep == "1"])
by_batch_subset$keep <- paste(by_batch_subset$keep, "QC_intrabatch_RSD_over30perc", sep = ";")        

##### Iterate over smaller data frame to update biomolecule keep values in sqlite db ######

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#dbExecute(con, "UPDATE biomolecules SET keep = '1'")                        

apply(by_batch_subset, 1, replaceValues)

## confirm 
biomolecule_df <- dbGetQuery(con, "SELECT *
          FROM biomolecules
          ")

dbDisconnect(con)

#### Establish a connection to the DB to obtain tier information #####

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

biomolecule_tier <- dbGetQuery(con, "SELECT *
                             FROM biomolecules
                             INNER JOIN metadata on metadata.biomolecule_id = biomolecules.biomolecule_id
                             WHERE omics_id = 3
                             AND (metadata_type = 'tier_mean') ", stringsAsFactors = F)

dbDisconnect(con)

#### create a subset of mean tier > 4 to flag #####

# create subset of features which need keep updated 
tier_subset <- data.frame( biomolecule_id = biomolecule_tier$biomolecule_id[biomolecule_tier$metadata_value > 4],keep = biomolecule_tier$keep[biomolecule_tier$metadata_value > 4] , stringsAsFactors = F)

# modify the keep string; if 1 turn to 0, then append "QC_2_intrabatch_RSD_over30perc"
tier_subset$keep[tier_subset$keep == "1"] <- sub("1", "0", tier_subset$keep[tier_subset$keep == "1"])
tier_subset$keep <- paste(tier_subset$keep, "mean_tier_greater_than_4", sep = ";")        

##### Update the database biomolecules table ######

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#dbExecute(con, "UPDATE biomolecules SET keep = '1'")                        

apply(tier_subset, 1, replaceValues)

## confirm 
biomolecule_df <- dbGetQuery(con, "SELECT *
          FROM biomolecules
          ")

dbDisconnect(con)
