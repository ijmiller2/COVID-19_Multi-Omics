###### 16_KAO_merging_CD3.1_results_with_Lipidex_output.R #######

## To help identifiy some of the unknown features we searched raw files
## with extra nodes in Compound Discoverer 3.1: MZ Cloud, chemspider, 
## mzvault, and formula.
## These results were further filtered down in excel to include only
## features (putative IDs) with matches to mz Cloud ("Full Match") or 
## mz Vault ("Full Match"). Next, because Compound discover output 
## collaspes the different adduct m/z into one compound. Another 
## 2 columns were added: m/z and adduct. If a compound had more than
## one adduct, a separate row was generated for each adduct with
## a copy of the identification results. This strategy, in theory, 
## should help with matching the unknowns from LipiDex output since
## Lipidex output report features m/z and RT. 
##
## This script is intendend to match unknown features by m/z and RT 
## between lipidex results and these modified CD3.1 output results. 
## The matching will be done by first rounding RT and m/z values and 
## then merging the 2 documents.


library(DBI)
library(RSQLite)

## load in lipid data from db

# connect 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# pull pvalues 
lipids <- dbGetQuery(con, "SELECT * FROM biomolecules
          WHERE omics_id = 2 AND keep = 1 
                  ")

# disconnect
dbDisconnect(con) 

## load in CD3.1 results file 

cd_results <- read.csv("P:/All_20200428_COVID_plasma_multiomics/Lipidomics/CD3_all_discovery_metabolomics_filtered.csv", stringsAsFactors = F)

names(cd_results)
cd_results <- cd_results[,1:28]

#### lipids standardized names contains mz and RT, need to round #### 

names(lipids)
lipids_unknowns <- lipids[grep("nknown", lipids$standardized_name),]

lipids_RT <- apply(lipids_unknowns, 1, function(x) unlist(strsplit(x[2], "RT_"))[2])
lipids_RT_round <- round(as.numeric(lipids_RT), digits = 2)

lipids_MZ <- apply(lipids_unknowns, 1, function(x) unlist(strsplit(unlist(strsplit(x[2], "mz_"))[2], "_"))[1])
lipids_MZ_round <- round(as.numeric(lipids_MZ), digits = 2)


## checking to see if any potential matches 
table(lipids_MZ_round %in% round(cd_results$m.z, digits =2 ))
length(lipids_MZ_round)
length(cd_results$m.z)      

table(lipids_RT_round %in% round(cd_results$RT..min., digits = 2))

cd_results$mz_RT <- paste(round(cd_results$m.z, digits = 2), round(cd_results$RT..min., digits = 2), sep ="_")

lipids_unknowns$mz_RT <- paste(lipids_MZ_round, lipids_RT_round, sep = "_")

##### merging two data sets#### 

merge_unknowns <- merge(lipids_unknowns, cd_results, by ="mz_RT")

write.csv(merge_unknowns, "data/Sup_table_2_merge_unknowns.csv")

##### Appending this information to the metadata table in DB #### 

## read current metdata table 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# pull pvalues 
metadata <- dbGetQuery(con, "SELECT * FROM metadata
                  ")

# disconnect
dbDisconnect(con) 

names(metadata)
df_metadata_append <- data.frame(metadata_id = NA, 
                                       biomolecule_id = merge_unknowns$biomolecule_id, 
                                       metadata_type = "Potential_annotation_through_secondary_db_searching", 
                                       metadata_value = merge_unknowns$Name)
                      
df_metadata_append$metadata_id <- seq(max(metadata$metadata_id)+1, length.out = nrow(df_metadata_append), by = 1)

## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

## write table to DB 

dbWriteTable(con, "metadata", df_metadata_append, append = T)

# check
metadata <- dbReadTable(con, "metadata")

# disconnect
dbDisconnect(con) 



