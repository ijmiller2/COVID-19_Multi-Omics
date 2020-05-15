##### X3_Updating_GC_metabolite_tier_in_DB.R #######

## Currently the biomolecule information in the sqlite db
## does not contain tier information, which would be 
## especially important to consider when filtering 
## metabolites to use for downstream analysis. 
## 
## The tier information is contained in the GC_quant 
## output file (.txt). This will be summarized and 
## added to the metadata table in the sqlite db
## 


library(DBI)
library(RSQLite)

##### Load in the GC data from P: Drive ##### 

df <- read.delim("P:/All_20200428_COVID_plasma_multiomics/GC_metabolomics/GC_Quant/QuantResults_202005111640419156_YH_KAO_annotated3.txt", stringsAsFactors = F)

##### Summarize tier information ##### 

tier_mean <- rowMeans(df[,grepl("tier", names(df))])

#### Establish a connection to the DB #####

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

##### Extracting the biomolecule identifier from sqlite db ##### 

biomolecule_info <- dbGetQuery(con, "SELECT biomolecule_id, standardized_name
           FROM biomolecules
           WHERE omics_id = 3")

dbDisconnect(con)

##### Confirming match of biomolecule names ###### 

all(biomolecule_info$standardized_name == df$Feature.ID)
## for some reason there is a duplicate
biomolecule_info <- biomolecule_info[1:length(df$Feature.ID),]

##### Linking biomolecule_id to tier information ##### 

df_tier <- data.frame(biomolecule_id = biomolecule_info$biomolecule_id, tier_mean = tier_mean)

##### R Function to insert values ##### 

insertValues <- function(x) dbExecute(con, paste("INSERT INTO metadata (biomolecule_id, metadata_type, metadata_value) 
          VALUES (",unlist(x[1]),", 'tier_mean', ", unlist(x[2]), ")" , sep = ""))

#### Establish a connection to the DB #####

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### Append data into db using apply function ##### 

apply(df_tier, 1, insertValues)

#### Confirm append worked #### 

x <- dbReadTable(con, "metadata")

dbDisconnect(con)
