
##### Code from: eda/KAO/01_KAO_Establishing_connection_to_db_extracting_timeStamp #####
##### Accessing time stamp of raw files #####

library(DBI)
library(RSQLite)

#### Establish connection to SQLite db #####

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### List tables in db #####

dbGetInfo(con)
dbListTables(con)
dbReadTable(con, "rawfiles")

#### Fetch raw file time stamp for GC files ####

timeStamp <- dbGetQuery(con, "SELECT * FROM rawfiles WHERE ome_id=3")

dbDisconnect(con)
