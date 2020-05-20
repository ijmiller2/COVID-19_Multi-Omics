#### 01_EAT_exploring_KAO_04_GCdata.R ##### 

## Updates to the database help with feature filtering of the 
## GC metabolites. See "X4_KAO_Updating_biomolecules_keep_column
## _for_GC_metabolites.R" and "03_KAO_Exploring_GC_feature_quality.R"
## 


## Access DB SQlite where GC data lives to use in pheatmap. 
## The code was taken from KAO 

library(DBI)
library(RSQLite)
library(pheatmap)
library(plyr)
library(RColorBrewer)

##### Establish a connection with DB #####

con <- dbConnect(SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

##### Pull data from DB #####

dbListTables(con)
# [1] "biomolecules"                  "deidentified_patient_metadata" "lipidomics_measurements"       "lipidomics_runs"              
# [5] "metabolomics_measurements"     "metabolomics_runs"             "metadata"                      "omes"                         
# [9] "patient_metadata"              "patient_samples"               "proteomics_measurements"       "proteomics_runs"              
# [13] "rawfiles"                      "sqlite_sequence"              

dbReadTable(con, "biomolecules")[0,]
dbReadTable(con, "metadata")[1:10,]
dbReadTable(con, "metabolomics_measurements")[0,]
dbReadTable(con, "metabolomics_runs")[0,]
dbReadTable(con, "biomolecules")[0,]

##### Knit into a long dataframe #####
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

##### Data formatting #####

# Reshape data into wide format with first two columns containing meta data
df_wide <- reshape(df, timevar = "biomolecule_id", v.names = "normalized_abundance",
                   idvar = "unique_identifier", direction = "wide" )

# Subset meta data and break up string by "_" 
meta <- as.data.frame(stringr::str_split_fixed(df_wide$unique_identifier,"_",5))
meta$V5 <- meta$V4 #overwrite split information of sample name
meta$V6 <- df_wide$batch #add batch info column
colnames(meta) <- c("Date", "Sample.Type", "Disease.State", "Patient.ID.unique","Patient.ID.group","Batch") #rename meta columns

# The patient id's are used more than once.
# To make all id's unique, the make.unique function will add numbers to each duplicated entry
meta$Patient.ID.unique <- make.unique(as.character(meta$Patient.ID.unique))

# Fix typo in Disease state. Update factors and eleves
meta$Disease.State[which(meta$Disease.State == "Conrol")] <- "Control"
meta$Disease.State <- factor(meta$Disease.State,levels = c("Control","COVID","NONCOVID"))

##### Prepare for heatmap ######

# Subset dataframe to include abunance values only
metabolomics <- df_wide[,c(3:ncol(df_wide))]
# rename row names to Patient unique Id. (non-duplicated)
rownames(metabolomics) <- meta$Patient.ID.unique
# rename columns to only the metabolite unique identifier
colnames(metabolomics) <- sub(".*abundance.","",colnames(metabolomics))

# Annotations for rows

row_annotations <- meta[,c(1,3,6)]
rownames(row_annotations) <- meta$Patient.ID.unique
row_annotations$Date <- as.Date(row_annotations$Date)

my_colour = list(
  Date = c(`20200427` = "#F4F4F9", `20200428` = "#E1FBFB", `20200429` = "#B8DBD9", `20200430` = "#86D6E0",
           `20200506` = "#3C99B2", `20200507` = "#368BA4", `20200508` = "#18475C" ),
  Disease.State  = c(`Control` = "#F7F7F7", `COVID` = "#D0D3CA", `NONCOVID` = "#A3A79D"),
  Batch = c(`1` = "#BCF8EC", `2` = "#AED9E0", `3` = "#9FA0C3", `4`  = "#8B687F", `5` = "#7B435B", `6` = "#5C3344", `7` = "#1C093D"))

##### Plot Heatmap #####
scaleRYG <- colorRampPalette(c("#3C99B2","#E8CB2E","#EF2D00"), space = "rgb")(20)

pheatmap(metabolomics,
         color = scaleRYG,
         annotation_colors = my_colour,
         annotation_row = row_annotations,
         cluster_rows = T,
         #scale = "column",
         main = "GC-MS Small Molecule COVID-19 HeatMap")

