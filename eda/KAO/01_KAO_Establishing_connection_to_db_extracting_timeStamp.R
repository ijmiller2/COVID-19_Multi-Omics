##### Accessing time stamp of raw files #####

setwd("eda/KAO")

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
