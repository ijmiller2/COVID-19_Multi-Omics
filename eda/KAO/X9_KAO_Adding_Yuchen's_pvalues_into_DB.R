##### X9_KAO_Adding_Yuchens_pvalues_into_DB.R #### 

## Yuchen performed analysis on HFD for each biomolecue 
## those data are found in Rdata files in regression folder

library(DBI)
library(RSQLite)

##### read in p-values that are currently in DB ####

## Establish a connection to the DB
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 

## read Rdata files

load("E:/COVID-19_Multi-Omics/data/regression/df_pvalues_Hospital_free_days_45.RData")

df_pvalues_Hospital_free_days_45$test <- paste("ANOVA_", df_pvalues_Hospital_free_days_45$test)
df_pvalues_Hospital_free_days_45$comparison <- "Hospital_free_days_45"
df_pvalues_Hospital_free_days_45$pvalue_id <- seq(nrow(pvalues)+1, length.out = nrow(df_pvalues_Hospital_free_days_45), by =1 )


## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

## write table to DB 

dbWriteTable(con, "pvalues", df_pvalues_Hospital_free_days_45, append = T)

# check
pvalues <- dbReadTable(con, "pvalues")

# disconnect
dbDisconnect(con) 
