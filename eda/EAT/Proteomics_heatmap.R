#### 01_EAT_Proteomics_ForUnsupervisedAnalysis.R ##### 

## Access DB SQlite where Proteomics data. 
## Append patient meta information

library(DBI)
library(RSQLite)
library(pheatmap)
library(plyr)
library(RColorBrewer)
library(DelayedMatrixStats)
library(randomcoloR)
library(cluster)

Deane_proteins <- read.csv("H:/Projects/COVID19/Proteomics/Files/131ProteinsofInterest_UniProtID.csv", 
                                                       header = TRUE, sep = ",", stringsAsFactors = FALSE)

##### Establish a connection with DB #####

con <- dbConnect(SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

##### Pull data from DB #####

dbListTables(con)
# [1] "biomolecules"                  "deidentified_patient_metadata" "lipidomics_measurements"       "lipidomics_runs"              
# [5] "metabolomics_measurements"     "metabolomics_runs"             "metadata"                      "omes"                         
# [9] "patient_metadata"              "patient_samples"               "proteomics_measurements"       "proteomics_runs"              
# [13] "rawfiles"                      "sqlite_sequence"               

# Extract tables from SQLite database. The square prakets provide the first row or column names of table. 
dbReadTable(con, "biomolecules")[0,]
dbReadTable(con, "metadata")[1:10,]
dbReadTable(con, "metabolomics_measurements")[0,]
dbReadTable(con, "metabolomics_runs")[0,]
dbReadTable(con, "biomolecules")[0,]
dbReadTable(con, "deidentified_patient_metadata")[0,]

dbReadTable(con, "proteomics_measurements")[0,]
#[1] measurement_id       replicate_id         biomolecule_id       raw_abundance        normalized_abundance raw_ibaq             normalized_ibaq 
dbReadTable(con, "proteomics_runs")[0,]
#[1] replicate_id      rawfile_id        unique_identifier

##### Knit tables into a long dataframe #####

df <- dbGetQuery(con, "SELECT proteomics_runs.unique_identifier, proteomics_measurements.normalized_abundance, proteomics_measurements.biomolecule_id, batch 
           FROM proteomics_measurements
           INNER JOIN proteomics_runs ON proteomics_runs.replicate_id = proteomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = proteomics_runs.rawfile_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = proteomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1 
           AND ome_id = 1 
           AND biomolecules.keep = '1'
           ")

# Select Sample_label, Gender, ICU_1, Mech_Ventilation from deidentified_patient_metadata
deidentified_patient_meta <- dbGetQuery(con, "SELECT Sample_label, Gender, ICU_1, Mech_Ventilation
                 FROM deidentified_patient_metadata
                 ")


proteins_meta <- dbGetQuery(con, "SELECT standardized_name, metadata.biomolecule_id, metadata_type, metadata_value
                               FROM metadata
                               INNER JOIN biomolecules ON biomolecules.biomolecule_id = metadata.biomolecule_id
                               WHERE biomolecules.omics_id = 1
                               AND biomolecules.keep = 1
                               ")
dbDisconnect(con)

##### Data formatting #####

# Reshape data into wide format with first two columns containing meta data
# Column names is "normalized_abundance.X" where X is the biomolecule.id
df_wide <- reshape(df, timevar = "biomolecule_id", v.names = "normalized_abundance",
                   idvar = "unique_identifier", direction = "wide" )

# Match rawfiles to match "20200514_ES_HC_afterBatch1" -> "20200514_ES_COON_HC_afterBatch1_single-shot"
HC_rows <- grep("HC",df_wide$unique_identifier)

for(i in 1:length(HC_rows)){
  row_df <- HC_rows[i]
  old <- df_wide$unique_identifier[row_df]
  new <- gsub('^(.{11})(.*)$', '\\1_COON\\2', old)
  new <- gsub('^(.{31})(.*)$', '\\1_single-shot\\2', new)
  
  df_wide$unique_identifier[row_df] <- new
}



##### Meta formatting #####

## meta contains the relevant information provided in the rawfile name (aka unique_identifier)

# Subset meta data and break up string by "_" 
meta <- as.data.frame(stringr::str_split_fixed(df_wide$unique_identifier,"_",7))

# Fix typo in Disease state. Update factors and eleves
meta$V4 <- sub("19","", meta$V4)
meta$V4 <- as.factor(sub("-","",meta$V4))
meta$V4 <- sub("control","pooled_plasma",meta$V4)

batch_row <- grep("single",meta$V5)

for(i in 1:length(batch_row)){
  row_meta <- batch_row[i]
  v5 <- meta$V5[row_meta]
  v6 <- meta$V6[row_meta]
  
  meta$V5[row_meta] <- v6
  meta$V6[row_meta] <- v5
}

meta$V5 <- sub("*.atch","", meta$V5)
meta$V5 <- as.factor(sub("after","", meta$V5))


# Pad single digit numbers with a zero ex. 2 -> '02'
samples <- vector()
for (i in 1:nrow(meta)) {
  if(nchar(as.character(meta$V7[i]),type = "char")==1){
    one_sample <- paste0("0",meta$V7[i])
    samples <- append(samples,one_sample)
  }else{
    one_sample <- as.character(meta$V7[i])
    samples <- append(samples,one_sample)
  }
  
}

#overwrite sample_number padded with a zero if single digit
meta$V7 <- samples 
meta$V8 <- paste0(meta$V4,"_",meta$V7)
# The patient id's are used more than once.
# To make all id's unique, the make.unique function will add numbers to each duplicated entry 
meta$V9 <- make.unique(as.character(meta$V8))

#add batch info column
meta <- meta[,c(1,4,5,7,8,9)]
#rename meta columns
colnames(meta) <- c("Rawfile.Date.Stamp","Disease.State", "Batch",
                    "Padded.Patient.Number", "Sample_Label", 
                    "Sample_Label.unique") 


##### Protein meta formattting #####

proteins_meta <- proteins_meta[which(proteins_meta$metadata_type == "gene_name"),]
# Isolate only the first protein in protein group for Majority.protein.IDs and Protein.IDs
majorityProtein <- lapply(strsplit(proteins_meta$standardized_name, ";"), '[', 1)
proteins_meta$majority.protein <- vapply(majorityProtein, paste, collapse = ", ", character(1L))

##### Prepare for heatmap ######

# Subset dataframe to include abunance values only
proteomics <- df_wide[,c(3:ncol(df_wide))]
# rename row names to Patient unique Id. (non-duplicated)
rownames(proteomics) <- meta$Sample_Label.unique
# rename columns to only the metabolite unique identifier
colnames(proteomics) <- sub(".*abundance.","",colnames(proteomics))


# Annotations for patients
patient_annotation <- meta[,c(2,3)]
patient_annotation$Disease.State <- as.factor(patient_annotation$Disease.State)
rownames(patient_annotation) <- meta$Sample_Label.unique
#patient_annotation$Date <- as.Date(patient_annotation$Date)

##### Plot Heatmap #####
scaleRYG <- colorRampPalette(c("#3C99B2","#E8CB2E","#EF2D00"), space = "rgb")(20)

pheatmap(t(proteomics),
         color = scaleRYG,
         #annotation_colors = my_colour,
         annotation_col = patient_annotation,
         cluster_rows = T,
         #scale = "column",
         show_colnames = F,
         show_rownames = F,
         main = "Proteomics COVID-19 HeatMap")

#### Remove Controls ####

## Create a heatmap without Controls and append more patient metadata

##### Patient Meta formatting ####

# Remove controls from meta 
meta_patientONLY <- meta[-grep("HC",meta$Sample_Label.unique),]
meta_patientONLY <- meta_patientONLY[-grep("pooled",meta_patientONLY$Sample_Label.unique),]

# What patient meta information is missing?
PatientSample_Missing <- meta_patientONLY[which(!meta_patientONLY$Sample_Label %in% deidentified_patient_meta$Sample_label),]
Deident_patientSample_Missing <- deidentified_patient_meta[which(!deidentified_patient_meta$Sample_label %in% meta_patientONLY$Sample_Label),]

#Missing samples 
# 121  NONCOVID_18      M     1                1
# 124  NONCOVID_21      M     0                0
# 126  NONCOVID_23      F     0                0
# 128  NONCOVID_25      M     0                0

# match sample labels in meta and deidentified_patient_meta
meta_patientONLY <- meta_patientONLY[which(meta_patientONLY$Sample_Label.unique %in% deidentified_patient_meta$Sample_label),]

# reorder deidentified_patient_meta data to match sample label in meta_patientONLY
deidentified_patient_meta <- deidentified_patient_meta[-which(rownames(deidentified_patient_meta) %in% rownames(Deident_patientSample_Missing)),]
deidentified_patient_meta <- deidentified_patient_meta[order(match(deidentified_patient_meta$Sample_label,meta_patientONLY$Sample_Label)),]

# Merge deidentified_patient_meta with meta_patientONLY by sample label
meta_patientONLY_merged <- merge(meta_patientONLY,deidentified_patient_meta, by.x = "Sample_Label", by.y = "Sample_label" )

##### Data without Controls #####

# Remove rows that are connected with Control Data
proteomics_noControls <- proteomics[which(rownames(proteomics) %in% meta_patientONLY_merged$Sample_Label),]

# For metaboanalyst
Disease.state <- sub("_.*", "",rownames(proteomics_noControls))
metaboanalyst_df <- cbind(Disease.state,proteomics_noControls)
metaboanalyst_df <- t(metaboanalyst_df)
rownames(metaboanalyst_df) <- c("Disease.state",proteins_meta$majority.protein)
#write.csv(metaboanalyst_df, file = "H:/Projects/COVID19/COVID_proteomics_onlyPatients_normalized.csv")

# Annotations for patients
patient_annotation <- meta_patientONLY_merged[,c(3,4,7,8,9)]
patient_annotation$Gender <- factor(patient_annotation$Gender, levels = c("F","M"))

# Substitute all 0 -> NO and 1 -> Yes in ICU_1 
patient_annotation$ICU_1 <- as.character(patient_annotation$ICU_1)
patient_annotation$ICU_1 <- sub("0","No",patient_annotation$ICU_1)
patient_annotation$ICU_1 <- sub("1","Yes",patient_annotation$ICU_1)
patient_annotation$ICU_1 <- factor(patient_annotation$ICU_1, levels = c("No","Yes"))

# Substitute all 1 -> NO and 1 -> Yes in ICU_1 
patient_annotation$Mech_Ventilation <- as.character(patient_annotation$Mech_Ventilation)
patient_annotation$Mech_Ventilation <- sub("0","No",patient_annotation$Mech_Ventilation)
patient_annotation$Mech_Ventilation <- sub("1","Yes",patient_annotation$Mech_Ventilation)
patient_annotation$Mech_Ventilation <- factor(patient_annotation$Mech_Ventilation, levels = c("No","Yes"))

rownames(patient_annotation) <- meta_patientONLY_merged$Sample_Label

my_colour = list(
  Disease.State  = c( `COVID` = "#D0D3CA", `NONCOVID` = "#A3A79D"),
  Batch = c(`1` = "#FDF0FE", `2` = "#E8C5E9", `3` = "#9FA0C3", `4`  = "#8B687F", `5` = "#7B435B", `6` = "#5C3344", `7` = "#1C093D"),
  Gender = c(`F` = "#C0B9DD", `M` = "#234947"),
  ICU_1 = c(`No` = "#D0D3CA", `Yes` = "#368BA4"),
  Mech_Ventilation = c(`No` = "#D0D3CA", `Yes` = "#463F3A"))

#scaleRYG <- colorRampPalette(c("#3C99B2","#E8CB2E","#EF2D00"), space = "rgb")(20)
#scaleRYG <- colorRampPalette(c("#3C99B2","#ffffff","#EF2D00"), space = "rgb")(20)
pheatmap(t(proteomics_noControls),
         color = scaleRYG,
         annotation_colors = my_colour,
         annotation_col = patient_annotation,
         #annotation_row = metabolite_annotation,
         cluster_cols = T,
         #scale = "column",
         show_colnames = F,
         show_rownames = F,
         #scale = "row",
         main = "Proteomics COVID-19 patients only HeatMap")


#### Coagulation Cascade ####

coagulation_proteins_meta <- proteins_meta[which(proteins_meta$majority.protein %in% Deane_proteins$UniProt.ID),]

coagulation_proteomics <- proteomics_noControls[,which(colnames(proteomics_noControls) %in% coagulation_proteins_meta$biomolecule_id)]
colnames(coagulation_proteomics) <- coagulation_proteins_meta$majority.protein

pheatmap(t(coagulation_proteomics),
         color = scaleRYG,
         annotation_colors = my_colour,
         annotation_col = patient_annotation,
         #annotation_row = metabolite_annotation,
         cluster_cols = F,
         #scale = "column",
         show_colnames = F,
         show_rownames = F,
         #scale = "row",
         main = "Coagulation Proteins of Interest\n COVID-19 patients only HeatMap")
