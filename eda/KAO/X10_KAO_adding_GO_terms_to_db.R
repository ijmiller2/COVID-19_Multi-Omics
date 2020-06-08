#### X10_KAO_adding_GO_terms_to_db.R #### 

## Anji provided a table of uniprot IDs to GO terms. This script links those 
## uniprot_ids to biomolecule_ids and adds these data into the db for future 
## GO enrichment analysis. 

library(DBI)
library(RSQLite)

## get biomolecule_id_standardized dame from db 

## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# pull biomolecule table and metadata table
biomolecules <- dbGetQuery(con, "SELECT * FROM biomolecules WHERE omics_id = 1")
metadata <- dbGetQuery(con, "SELECT * FROM metadata")

# disconnect
dbDisconnect(con) 

#### read in Anji's table ##### 

GO_terms <- read.csv("P:/All_20200428_COVID_plasma_multiomics/Proteomics/uniprot_GO_BioProcess_CellularComp_MolecularFun.csv", stringsAsFactors = F)

dim(GO_terms) #513
names(GO_terms)

GO_terms$Gene.ontology..biological.process.<- sub("[", "", GO_terms$Gene.ontology..biological.process. , fixed =T )
GO_terms$Gene.ontology..biological.process.<- sub("]", "", GO_terms$Gene.ontology..biological.process. , fixed =T )

GO_terms$Gene.ontology..cellular.component. <- sub("[", "", GO_terms$Gene.ontology..cellular.component. , fixed =T )
GO_terms$Gene.ontology..cellular.component.<- sub("]", "", GO_terms$Gene.ontology..cellular.component. , fixed =T )

GO_terms$Gene.ontology..molecular.function. <- sub("[", "", GO_terms$Gene.ontology..molecular.function. , fixed =T )
GO_terms$Gene.ontology..molecular.function. <- sub("]", "", GO_terms$Gene.ontology..molecular.function. , fixed =T )

GO_terms$yourlist <- apply(GO_terms, 1, function(x) strsplit(x[8], ";")[[1]][1])

#### extracting out first uniprot_id from standardized_name to match with Anji's table #### 

first<- apply(biomolecules, 1, function(x) strsplit(x[2], ";")[[1]][1])

biomolecules$first <-first 
any(duplicated(biomolecules$first ))

#### Merge Anji's table with biomoleucles ###### 

merge_df <- merge(GO_terms, biomolecules, by.x = "yourlist", by.y = "first", all.x =T)

#### Get info from metadata table in order to append GO term information ##### 

names(metadata) #46189

df_metadata_append <- rbind(data.frame(metadata_id = NA, 
                                       biomolecule_id = merge_df$biomolecule_id, 
                                       metadata_type = "GO_biological_process", 
                                       metadata_value = merge_df$Gene.ontology..biological.process.),
                            data.frame(metadata_id = NA,
                                       biomolecule_id = merge_df$biomolecule_id,
                                       metadata_type = "GO_cellular_component",
                                       metadata_value = merge_df$Gene.ontology..cellular.component.),
                            data.frame(metadata_id = NA,
                                       biomolecule_id = merge_df$biomolecule_id,
                                       metadata_type = "GO_molecular_function",
                                       metadata_value = merge_df$Gene.ontology..molecular.function.))

df_metadata_append$metadata_id <- seq(max(metadata$metadata_id)+1, length.out = nrow(df_metadata_append), by = 1)

#### Add GO terms to metadata #### 
## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### write table to DB 

dbWriteTable(con, "metadata", df_metadata_append, append = T)

# check
metadata <- dbReadTable(con, "metadata")

# disconnect
dbDisconnect(con) 

