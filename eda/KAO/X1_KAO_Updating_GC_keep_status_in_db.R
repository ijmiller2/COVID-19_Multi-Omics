##### MODIFYING DB 2020-05-13 #####

## edit by KAO
## GC files in Batch 1 were run with 1:10 split and 1:25 split,
## but 1:25 split resulted in better peak shapes and was considered
## more optimal for the future runs. 

## rawfiles table in sqlite db has column 'keep' and is
## set as 1 (TRUE) for all files. This script turns Batch 1 1:10 
## split files to keep = 0 (FALSE).

library(DBI)
library(RSQLite)


#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")


##### Test selection of batch 1 raw files with 1:10 split ##### 
dbGetQuery(con, "SELECT keep 
          FROM rawfiles
          WHERE rawfile_name LIKE '%20200427_BP%' 
          AND rawfile_name NOT LIKE '%25split%'
          ")

#### Set keep to 0 for these selected entries ####
dbExecute(con, "UPDATE rawfiles
          SET keep = 0
          WHERE rawfile_name LIKE '%20200427_BP%' 
          AND rawfile_name NOT LIKE '%25split%'
          ")

#### Confirm that keep is now 0 for these files ####
dbGetQuery(con, "SELECT *
          FROM rawfiles
          WHERE rawfile_name LIKE '%20200427_BP%' 
          AND rawfile_name NOT LIKE '%25split%'
          ")

dbGetQuery(con, "SELECT *
          FROM rawfiles
          WHERE keep = 0
          ")

dbDisconnect(con)
