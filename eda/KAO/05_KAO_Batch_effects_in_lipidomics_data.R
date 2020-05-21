##### 05_KAO_Lipidomics_batch_effects.R 

## This script is to see if the lipidomics data need to be corrected for batch
## effects. 

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
dbReadTable(con, "lipidomics_measurements")[0,]
dbReadTable(con, "lipidomics_runs")
dbReadTable(con, "biomolecules")[0,]


##### Creating df of lipid raw abundance ###### 
lipids <- dbGetQuery(con, "SELECT unique_identifier, timestamp, raw_abundance, biomolecule_id, batch
           FROM lipidomics_measurements
           INNER JOIN lipidomics_runs ON lipidomics_runs.replicate_id = lipidomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = lipidomics_runs.rawfile_id
           ")

lipid_ids <- dbGetQuery(con, "SELECT * 
                        FROM biomolecules
                        WHERE omics_id = 2")

dbDisconnect(con)

#### Looking for potential batch effect in lipids data #### 

## PLOT1: boxplot of batch intensity for all samples 
boxplot(log2(lipids$raw_abundance)~lipids$batch, 
        col = 1:7, 
        ylab = "Lipid feature intensity log2()",
        xlab = "Batch",
        main = "Lipid feature intensity by batch")


controls <- grepl("Control", lipids$unique_identifier) & !grepl("Control_half", lipids$unique_identifier) & !grepl("Control_quarter", lipids$unique_identifier)

## PLOT2: boxplot of control lipid intensity across batches
boxplot(log2(lipids$raw_abundance[controls])~lipids$unique_identifier[controls],
        col = rep(1:7, each = 4),
        las = 2,
        main = "Lipid feature intensity of controls",
        ylab = "Lipid feature intensity log2() ")

##### Priciple component of lipids data ##### 

## convert to wide 
pre <- reshape(lipids, timevar = "biomolecule_id", v.names = "raw_abundance",
               idvar = "unique_identifier", direction = "wide" )

# selecting only data which is either controls or samples
pre <- pre[!(grepl("resolvin|Method", pre$unique_identifier)|grepl("half|quarter", pre$unique_identifier)), ]

pre[!grepl("Control", pre$unique_identifier),]$unique_identifier


# performing principle component analysis on raw abundance data 
pca_pre <- prcomp(log2(pre[,-c(1:3)]))

## PLOT3: PCA of raw lipid features #### 
plot(pca_pre$x, col = pre$batch, pch = 19, main = "PCA lipidomics - color by batch")

# controls appear to have batch effect, exploring the loadings for control-only pca
pca_pre_controls <- prcomp(log2(pre[grepl("Control", pre$unique_identifier),-c(1:3)]))

## PLOT4; PCA of controls only
plot(pca_pre_controls)
plot(pca_pre_controls$x, col= pre$batch[grepl("Control", pre$unique_identifier)], pch =19,
     main = "PCA lipids - raw abundance controls only")

## PLOT5 loadings in PC-1 controls only 
barplot(pca_pre_controls$rotation[,1][order(pca_pre_controls$rotation[,1])])

pca_pre_controls$rotation[,1]        
dim(pre)
dim(lipid_ids)

#### Matching lipids to lipid_ids #####
range(unique(lipids$biomolecule_id)) # 225-7480
range(lipid_ids$biomolecule_id) # 216-7480

length(unique(lipids$biomolecule_id)) #6814
length(unique(lipid_ids$biomolecule_id)) #7265
length(unique(lipid_ids$standardized_name)) #6841 
#means duplicate names, but even that doesn't make numbers match. 

table(lipid_ids$keep) # all are keep right now. 

length(grep("Duplicate", lipid_ids$standardized_name)) #1384

#### Figure out why not matching #####

## load in original lipid file. 
lipids_original <- read.csv("P:/All_20200428_COVID_plasma_multiomics/Lipidomics/Lipidomics_quant_results/Final_Results.csv", stringsAsFactors = F)

## load in the results file that Dain used
lipids_original_copy <- read.csv("P:/All_20200428_COVID_plasma_multiomics/Lipidomics/Lipidomics_quant_results/Final_Results (copied because file was locked).csv", stringsAsFactors = F)

## check dimensions
dim(lipids_original) #7235 191  
dim(lipids_original_copy) #7235 191 # these two file appear to be the same, but do not match other dimentions
dim(lipid_ids) #7265 5
length(unique(lipids$biomolecule_id)) #6814

### trying to match lipid_ids to lipid original names to figure out what's going on

lipid_ids$standardized_name[1:10] # first 10 match
lipids_original[1:10, 1:5]

lipid_ids$standardized_name[100:110] # 100- 110 do not match! 
lipids_original[100:110, 1:5]

# no finding between 10 and 100 where things don't match
lipid_ids$standardized_name[10:20] # 10-20 match
lipids_original[10:20, 1:5]

lipid_ids$standardized_name[20:30] # 20-30 match
lipids_original[20:30, 1:5]

lipid_ids$standardized_name[30:40] # 30 matches, but 31 does not match! Looks like it starts from beginning
lipids_original[30:40, 1:5]

## the dimension difference between original table and 
## sqlite db table is exactly 30, suggesting these first
## 30 are duplicates. 

lipid_ids2 <- lipid_ids[-c(1:30), ]
length(unique(lipid_ids2$standardized_name)) # 6814 
length(unique(lipids$biomolecule_id)) # 6814
# and now remaining duplicates are removed from the biomolecues table 
# to account for the descrpancy. Will need to update the db. 
