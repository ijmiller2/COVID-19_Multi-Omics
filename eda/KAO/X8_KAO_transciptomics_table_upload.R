####### 08_KAO_Getting transcriptomics data ready for DB ######### 

## transcriptomics data is now in P: Drive, these are normalized
## by median counts and are appropriate for between samples comparisons
## for the db, need to add rawfiles index: 

library(DBI)
library(RSQLite)
###### Need to reference biomolecule_id #########

## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# check
biomolecules <- dbReadTable(con, "biomolecules")

# disconnect
dbDisconnect(con) 

###### Read in transcriptomics table and organize into datatables for DB ####### 

df <- read.delim("P:/All_20200428_COVID_plasma_multiomics/Transcriptomics/genes.ec.no_hg.norm.tsv", stringsAsFactors = F)
df_transformed <- read.delim("P:/All_20200428_COVID_plasma_multiomics/Transcriptomics/genes.l2ec.no_hg.norm.tsv", stringsAsFactors = F)

dim(df)
names(df)
###### creating rawfiles table from names of samples ######

rawfiles <- data.frame(rawfile_id = seq(771, length.out = ncol(df)-1, by = 1),
                       timestamp = -1,
                       rawfile_name = names(df)[-1],
                       sample_id = c(1:53,55:80, 82:87, 89:128),
                       run_type = "Sample",
                       keep = 1, 
                       batch = 1,
                       ome_id = 5, stringsAsFactors = F
                       )

transcriptomics_run <- data.frame(replicate_id = c(1:ncol(df[,-1])), 
                                  rawfile_id = rawfiles$rawfile_id,
                                  unique_identifier = rawfiles$rawfile_name, stringsAsFactors = F)


biomolecules_append <- data.frame(biomolecule_id = seq(max(biomolecules$biomolecule_id + 1), length.out = nrow(df), by = 1),
                                  standardized_name = df$symbol, 
                                  omics_id = 5,
                                  keep = 1)

row.names(df) <- biomolecules_append$biomolecule_id


###### converting original to long format ######

df_long <- reshape(df, idvar = "standardized_name", 
                                ids = df$symbol, 
                                times = names(df)[-c(1)], 
                                timevar = "unique_identifier",
                                varying = list(names(df)[-c(1)]),
                                direction = "long")
head(df_long)

df_long_transformed <- reshape(df_transformed, idvar = "standardized_name", 
                               ids = df$symbol, 
                               times = names(df)[-c(1)], 
                               timevar = "unique_identifier",
                               varying = list(names(df)[-c(1)]),
                               direction = "long")

transcriptomics_measurements <- data.frame(measurement_id = 1:nrow(df_long), 
                                           replicate_id = rep(transcriptomics_run$replicate_id, each = nrow(biomolecules_append)),
                                           biomolecule_id = rep(biomolecules_append$biomolecule_id, nrow(transcriptomics_run)), 
                                           raw_abundance = -1, 
                                           normalized_abundance = df_long_transformed$C001,
                                           normalized_counts = df_long$C001)

head(transcriptomics_measurements)

#### changing keep status for transcripts 

hist(rowSums(df[,-1] > 10), main = "number of samples\nwhere value is > 0")

biomolecules_append$keep[rowSums(df[,-1]> 10) < 125*0.25] <- "0;75perc_samples_have_less_than_10_counts"  
biomolecules_append$keep

hist(transcriptomics_measurements$normalized_counts)

##### Adding these transcriptomics tables to the db ######

## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

## Append rawfiles 

dbWriteTable(con, "rawfiles", rawfiles, append = T)

## Append biomolecules table 

dbWriteTable(con, "biomolecules", biomolecules_append, append = T)

## Add transcriptomics_runs

dbWriteTable(con, "transcriptomics_runs", transcriptomics_run)

## Add transcriptomics_measurements 

dbWriteTable(con, "transcriptomics_measurements", transcriptomics_measurements)

dbDisconnect(con)




