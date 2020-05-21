#### 01_EAT_exploring_KAO_04_GCdata.R ##### 

## Access DB SQlite where GC data lives and use in pheatmap. 
## Append patient meta information
## The RSQLite code was taken from KAO.

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

# Extract tables from SQLite database. The square prakets provide the first row or column names of table. 
dbReadTable(con, "biomolecules")[0,]
dbReadTable(con, "metadata")[1:10,]
dbReadTable(con, "metabolomics_measurements")[0,]
dbReadTable(con, "metabolomics_runs")[0,]
dbReadTable(con, "biomolecules")[0,]
dbReadTable(con, "deidentified_patient_metadata")[0,]

##### Knit into a long dataframe #####

# Select unique_identifier from metabolomics_runs
# Select normalized_abundance from from metabolomics_measurements
# Select biomolecule.id from metabolomic_measurements
df<- dbGetQuery(con, "SELECT unique_identifier, normalized_abundance, metabolomics_measurements.biomolecule_id, batch 
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1 
           AND ome_id = 3 
           AND biomolecules.keep = '1'
           ")

deidentified_patient_meta <- dbGetQuery(con, "SELECT Sample_label, Gender, ICU_1, Mech_Ventilation
                 FROM deidentified_patient_metadata
                 ")
dbDisconnect(con)


##### Data formatting #####

# Reshape data into wide format with first two columns containing meta data
df_wide <- reshape(df, timevar = "biomolecule_id", v.names = "normalized_abundance",
                   idvar = "unique_identifier", direction = "wide" )


# Subset meta data and break up string by "_" 
meta <- as.data.frame(stringr::str_split_fixed(df_wide$unique_identifier,"_",5))
# Fix typo in Disease state. Update factors and eleves
meta$V3[which(meta$V3 == "Conrol")] <- "Control"
meta$V3 <- factor(meta$V3,levels = c("Control","COVID","NONCOVID"))

samples <- vector()
for (i in 1:nrow(meta)) {
  if(nchar(as.character(meta$V4[i]),type = "char")==1){
    one_sample <- paste0("0",meta$V4[i])
    samples <- append(samples,one_sample)
  }else{
    one_sample <- as.character(meta$V4[i])
    samples <- append(samples,one_sample)
  }
  
}


meta$V5 <- samples #overwrite split1:25 with sample_number padded with a zero if single digit
meta$V6 <- paste0(meta$V3,"_",meta$V5) 


meta$V7 <- df_wide$batch #add batch info column
colnames(meta) <- c("Rawfile.Date.Stamp", "Sample.Type", "Disease.State", "Patient.Number","Padded.Patient.Number", "Sample_Label", "Batch") #rename meta columns

# The patient id's are used more than once.
# To make all id's unique, the make.unique function will add numbers to each duplicated entry
meta$Sample_Lable.unique <- make.unique(as.character(meta$Sample_Label))


##### Prepare for heatmap ######

# Subset dataframe to include abunance values only
metabolomics <- df_wide[,c(3:ncol(df_wide))]
# rename row names to Patient unique Id. (non-duplicated)
rownames(metabolomics) <- meta$Sample_Lable.unique
# rename columns to only the metabolite unique identifier
colnames(metabolomics) <- sub(".*abundance.","",colnames(metabolomics))

# Annotations for rows
row_annotations <- meta[,c(1,3,7)]
row_annotations$Batch <- factor(row_annotations$Batch, levels = seq(1:7))
rownames(row_annotations) <- meta$meta$Sample_Lable.unique
#row_annotations$Date <- as.Date(row_annotations$Date)

my_colour = list(
  Rawfile.Date.Stamp = c(`20200427` = "#F4F4F9", `20200428` = "#E1FBFB", `20200429` = "#B8DBD9", `20200430` = "#86D6E0",
           `20200506` = "#3C99B2", `20200507` = "#368BA4", `20200508` = "#18475C" ),
  Disease.State  = c(`Control` = "#F7F7F7", `COVID` = "#D0D3CA", `NONCOVID` = "#A3A79D"),
  Batch = c(`1` = "#BCF8EC", `2` = "#AED9E0", `3` = "#9FA0C3", `4`  = "#8B687F", `5` = "#7B435B", `6` = "#5C3344", `7` = "#1C093D"))


##### Plot Heatmap #####
scaleRYG <- colorRampPalette(c("#3C99B2","#E8CB2E","#EF2D00"), space = "rgb")(20)

pheatmap(t(metabolomics),
         color = scaleRYG,
         annotation_colors = my_colour,
         annotation_col = row_annotations,
         cluster_rows = T,
         #scale = "column",
         show_colnames = F,
         show_rownames = F,
         main = "GC-MS Small Molecule COVID-19 HeatMap")

#### Remove Controls ####
meta_patientONLY <- meta[-grep("Control",meta$Sample_Lable.unique),]

PatientSample_Missing <- meta_patientONLY[which(!meta_patientONLY$Sample_Label %in% deidentified_patient_meta$Sample_label),]

meta_patientONLY <- meta_patientONLY[which(meta_patientONLY$Sample_Lable.unique %in% deidentified_patient_meta$Sample_label),]

deidentified_patient_meta <- deidentified_patient_meta[order(match(deidentified_patient_meta$Sample_label,meta_patientONLY$Sample_Label)),]

meta_patientONLY_merged <- merge(meta_patientONLY,deidentified_patient_meta, by.x = "Sample_Label", by.y = "Sample_label" )

# Remove rows that are connected with Control Datat
metabolomics_noControls <- metabolomics[which(rownames(metabolomics) %in% meta_patientONLY_merged$Sample_Label),]


# Annotations for rows
patient_annotation <- meta_patientONLY_merged[,c(2,4,7,9,10,11)]
patient_annotation$Batch <- factor(patient_annotation$Batch, levels = seq(1:7))
patient_annotation$Gender <- factor(patient_annotation$Gender, levels = c("F","M"))
patient_annotation$ICU_1 <- factor(patient_annotation$ICU_1, levels = c("0","1"))
patient_annotation$Mech_Ventilation <- factor(patient_annotation$Mech_Ventilation, levels = c("0","1"))

rownames(patient_annotation) <- meta_patientONLY_merged$Sample_Label
#patient_annotation$Date <- as.Date(patient_annotation$Date)

my_colour = list(
  Rawfile.Date.Stamp = c(`20200427` = "#F4F4F9", `20200428` = "#E1FBFB", `20200429` = "#B8DBD9", `20200430` = "#86D6E0",
                         `20200506` = "#3C99B2", `20200507` = "#368BA4", `20200508` = "#18475C" ),
  Disease.State  = c( `COVID` = "#D0D3CA", `NONCOVID` = "#A3A79D"),
  Batch = c(`1` = "#FDF0FE", `2` = "#E8C5E9", `3` = "#9FA0C3", `4`  = "#8B687F", `5` = "#7B435B", `6` = "#5C3344", `7` = "#1C093D"),
  Gender = c(`F` = "#FB3640", `M` = "#1D3461"))
  #ICU_1 = c(`0` = "#E6F3F6 ", `1` = "#63B0C6"),
  #Mech_Ventilation = c(`0` = "#D9FF29", `1` = "#1B2D2A"))

scaleRYG <- colorRampPalette(c("#3C99B2","#E8CB2E","#EF2D00"), space = "rgb")(20)
pheatmap(t(metabolomics_noControls),
         color = scaleRYG,
         annotation_colors = my_colour,
         annotation_col = patient_annotation,
         #cluster_rows = T,
         #scale = "column",
         show_colnames = F,
         show_rownames = F,
         main = "GC-MS Small Molecule COVID-19 HeatMap")
