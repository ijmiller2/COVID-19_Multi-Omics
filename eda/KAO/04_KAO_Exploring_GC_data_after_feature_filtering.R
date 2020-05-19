#### 04_KAO_exploring_GC_data_after_feature_filtering.R ##### 

## Updates to the database help with feature filtering of the 
## GC metabolites. See "X4_KAO_Updating_biomolecules_keep_column
## _for_GC_metabolites.R" and "03_KAO_Exploring_GC_feature_quality.R"
## 
## I would like to see how this feature filtering changes the 
## pca plots and what classes of molecules remain after 
## filtering. 


library(DBI)
library(RSQLite)
colors <- read.csv("../../reference/color_palette.txt",  stringsAsFactors = F)[,2]
palette(colors)

#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### Pull data from DB ####

dbListTables(con)
dbReadTable(con, "biomolecules")[0,]
dbReadTable(con, "metadata")[1:10,]
dbReadTable(con, "metabolomics_measurements")[0,]
dbReadTable(con, "metabolomics_runs")[0,]
dbReadTable(con, "biomolecules")[0,]


##### 
df<- dbGetQuery(con, "SELECT unique_identifier, normalized_abundance, metabolomics_measurements.biomolecule_id, batch 
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1 
           AND ome_id = 3 
           AND biomolecules.keep = '1'
           ")

dbDisconnect(con)


#### Creating a wide-format data frame to facilitate PCA #####

df_wide <- reshape(df, timevar = "biomolecule_id", v.names = "normalized_abundance",
               idvar = "unique_identifier", direction = "wide" )

df_wide_exprs <- df_wide[,-c(1:2)]

##### Priciple component analysis #######
pca1 <- prcomp(df_wide_exprs)

plot(pca1)
plot(pca1$x, col = df_wide$batch, pch =19)
