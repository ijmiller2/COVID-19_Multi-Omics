####### X11_Adding_GO_terms_for_transcripts_into_db.R ##### 

## I used uniprot to pull GO information for the transcripts.
## They were stored as a file but need to be loaded in an linked to 
## biomolecule_ids. 

library(DBI)
library(RSQLite)


#### reading in biomolecules table from db for transcripts #####

# connect 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# pull pvalues 
biomolecules <- dbGetQuery(con, "SELECT * FROM biomolecules 
          WHERE omics_id = 5 AND keep = 1
                  ")

metadata <- dbGetQuery(con, "SELECT * FROM metadata")


# disconnect
dbDisconnect(con) 

#### reading in the the data ######

go_genes <- read.delim("data/uniprot-genelist.tab")

head(go_genes)
names(go_genes)

go_gene_subset <- data.frame("standardized_name" = go_genes$yourlist.M20200625A94466D2655679D1FD8953E075198DA808E5F76, 
                             "GO_biological_process" = go_genes$Gene.ontology..biological.process.,
                             "GO_cellular_component" = go_genes$Gene.ontology..cellular.component.,
                             "GO_molecular_function" = go_genes$Gene.ontology..molecular.function.)

##### mering uniprot output with biomolecule_id #####
merge_go_biomolecules <- merge(biomolecules, go_gene_subset, by = "standardized_name")


######

df_metadata_append <- data.frame(metadata_id = NA, 
                                       biomolecule_id = merge_go_biomolecules$biomolecule_id, 
                                       metadata_type = "GO_biological_process", 
                                       metadata_value = merge_go_biomolecules$GO_biological_process)
                      
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

