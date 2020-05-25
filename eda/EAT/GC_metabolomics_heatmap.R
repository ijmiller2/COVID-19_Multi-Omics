#### 01_EAT_exploring_KAO_04_GCdata.R ##### 

## Access DB SQlite where GC data lives and use in pheatmap. 
## Append patient meta information
## The RSQLite code was taken from KAO.

BiocManager::install("DelayedMatrixStats")

library(DBI)
library(RSQLite)
library(pheatmap)
library(plyr)
library(RColorBrewer)
library(DelayedMatrixStats)
library(randomcoloR)

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

##### Knit tables into a long dataframe #####

## df is a dataframe of normalized metabolite abundances, the rawfile where metabolite was extracted from, and batch info
## deiidentified_patient_meta is a dataframe that contains metat data for the patients
## metabolite_class is a dataframe of the metabolite standarized names and metadata associated to metabolite class

# Select unique_identifier from metabolomics_runs
# Select normalized_abundance from from metabolomics_measurements
# Select biomolecule.id from metabolomic_measurements
# Select batch from rawfiles
df<- dbGetQuery(con, "SELECT unique_identifier, normalized_abundance, metabolomics_measurements.biomolecule_id, batch 
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1 
           AND ome_id = 3 
           AND biomolecules.keep = '1'
           ")

# Select Sample_label, Gender, ICU_1, Mech_Ventilation from deidentified_patient_metadata
deidentified_patient_meta <- dbGetQuery(con, "SELECT Sample_label, Gender, ICU_1, Mech_Ventilation
                 FROM deidentified_patient_metadata
                 ")


# Select Metabolite standardized_name in biomolecules
# Select biomolecule_id in metadatat
# Select metadata_type in metadata
# Select metadata_value in metadata
metabolite_class <- dbGetQuery(con, "SELECT standardized_name, metadata.biomolecule_id, metadata_type, metadata_value
                               FROM metadata
                               INNER JOIN biomolecules ON biomolecules.biomolecule_id = metadata.biomolecule_id
                               WHERE biomolecules.omics_id = 3
                               AND biomolecules.keep = 1
                               AND metadata.metadata_type = 'Class'
                               ")

dbDisconnect(con)



##### Data formatting #####

# Reshape data into wide format with first two columns containing meta data
# Column names is "normalized_abundance.X" where X is the biomolecule.id
df_wide <- reshape(df, timevar = "biomolecule_id", v.names = "normalized_abundance",
                   idvar = "unique_identifier", direction = "wide" )


##### Meta formatting #####

## meta contains the relevant information provided in the rawfile name (aka unique_identifier)

# Subset meta data and break up string by "_" 
meta <- as.data.frame(stringr::str_split_fixed(df_wide$unique_identifier,"_",5))
# Fix typo in Disease state. Update factors and eleves
meta$V3[which(meta$V3 == "Conrol")] <- "Control"
meta$V3 <- factor(meta$V3,levels = c("Control","COVID","NONCOVID"))

# Pad single digit numbers with a zero ex. 2 -> '02'
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

#overwrite split1:25 with sample_number padded with a zero if single digit
meta$V5 <- samples 
meta$V6 <- paste0(meta$V3,"_",meta$V5) 

#add batch info column
meta$V7 <- df_wide$batch
#rename meta columns
colnames(meta) <- c("Rawfile.Date.Stamp", "Sample.Type", 
                    "Disease.State", "Patient.Number",
                    "Padded.Patient.Number", "Sample_Label", 
                    "Batch") 

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

# Annotations for patients
patient_annotation <- meta[,c(1,3,7)]
patient_annotation$Batch <- factor(patient_annotation$Batch, levels = seq(1:7))
rownames(patient_annotation) <- meta$meta$Sample_Lable.unique
#patient_annotation$Date <- as.Date(patient_annotation$Date)

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
         annotation_col = patient_annotation,
         cluster_rows = T,
         #scale = "column",
         show_colnames = F,
         show_rownames = F,
         main = "GC-MS Small Molecule COVID-19 HeatMap")

#### Remove Controls ####

## Create a heatmap without Controls and append more patient metadata

##### Patient Meta formatting ####

# Remove controls from meta 
meta_patientONLY <- meta[-grep("Control",meta$Sample_Lable.unique),]

# What patient meta information is missing?
PatientSample_Missing <- meta_patientONLY[which(!meta_patientONLY$Sample_Label %in% deidentified_patient_meta$Sample_label),]

# match sample labels in meta and deidentified_patient_meta
meta_patientONLY <- meta_patientONLY[which(meta_patientONLY$Sample_Lable.unique %in% deidentified_patient_meta$Sample_label),]

# reorder deidentified_patient_meta data to match sample label in meta_patientONLY
deidentified_patient_meta <- deidentified_patient_meta[order(match(deidentified_patient_meta$Sample_label,meta_patientONLY$Sample_Label)),]

# Merge deidentified_patient_meta with meta_patientONLY by sample label
meta_patientONLY_merged <- merge(meta_patientONLY,deidentified_patient_meta, by.x = "Sample_Label", by.y = "Sample_label" )


##### Data without Controls #####

# Remove rows that are connected with Control Data
metabolomics_noControls <- metabolomics[which(rownames(metabolomics) %in% meta_patientONLY_merged$Sample_Label),]

# For metaboanalyst
Disease.state <- sub("_.*", "",rownames(metabolomics_noControls))
metaboanalyst_df <- cbind(Disease.state,metabolomics_noControls)
metaboanalyst_df <- t(metaboanalyst_df)
rownames(metaboanalyst_df) <- c("Disease.state",metabolite_class$standardized_name)
#write.csv(metaboanalyst_df, file = "H:/Projects/COVID19/SmallMolecule/Files/metaboanalyst_GCMetabolomics_normalized.csv")

# Annotations for patients
#patient_annotation <- meta_patientONLY_merged[,c(2,4,7,9,10,11)]
patient_annotation <- meta_patientONLY_merged[,c(4,9,10,11)]
#patient_annotation$Batch <- factor(patient_annotation$Batch, levels = seq(1:7))
patient_annotation$Gender <- factor(patient_annotation$Gender, levels = c("F","M"))
patient_annotation$ICU_1 <- factor(patient_annotation$ICU_1, levels = c("0","1"))
patient_annotation$Mech_Ventilation <- factor(as.character(patient_annotation$Mech_Ventilation), levels = c("0","1"))

rownames(patient_annotation) <- meta_patientONLY_merged$Sample_Label
#patient_annotation$Date <- as.Date(patient_annotation$Date)

# Annotations for molecules
metabolite_annotation  <- metabolite_class[which(metabolite_class$biomolecule_id %in% colnames(metabolomics)),]
rownames(metabolite_annotation) <- metabolite_annotation$biomolecule_id

metabolite_class_extended <- as.data.frame(stringr::str_split_fixed(metabolite_class$metadata_value,";",2))
metabolite_class_extended$V1 <- as.character(metabolite_class_extended$V1)
metabolite_class_extended$V2 <- trimws(as.character(metabolite_class_extended$V2))

simplified_class <- vector()
for (i in 1:nrow(metabolite_class_extended)) {
  if(nchar(as.character(metabolite_class_extended$V2[i]), type = "char") == 0){
    metabolite_class_extended$V2[i] <- metabolite_class_extended$V1[i]
  }
  
  if(length(grep("Amino",metabolite_class_extended$V2[i])) != 0){
    simplified_class <- append(simplified_class,"Amino acids")
  }else{
    simplified_class <- append(simplified_class,metabolite_class_extended$V2[i])
  }
}

metabolite_class_extended$V3 <- simplified_class


metabolite_annotation <- as.data.frame(metabolite_class_extended$V3)
colnames(metabolite_annotation) <- "Metabolite.Class"

metabolite_annotation$Metabolite.Class <- as.factor(metabolite_annotation$Metabolite.Class)




my_colour = list(
  #Rawfile.Date.Stamp = c(`20200427` = "#F4F4F9", `20200428` = "#E1FBFB", `20200429` = "#B8DBD9", `20200430` = "#86D6E0",
                         #`20200506` = "#3C99B2", `20200507` = "#368BA4", `20200508` = "#18475C" ),
  Disease.State  = c( `COVID` = "#D0D3CA", `NONCOVID` = "#A3A79D"),
  #Batch = c(`1` = "#FDF0FE", `2` = "#E8C5E9", `3` = "#9FA0C3", `4`  = "#8B687F", `5` = "#7B435B", `6` = "#5C3344", `7` = "#1C093D"),
  Gender = c(`F` = "#FB3640", `M` = "#1D3461"))
  #ICU_1 = c(`0` = "#E1FBFB ", `1` = "#368BA4"))
  #Mech_Ventilation = c(`0` = "#D9FF29", `1` = "#1B2D2A"))

#scaleRYG <- colorRampPalette(c("#3C99B2","#E8CB2E","#EF2D00"), space = "rgb")(20)
#scaleRYG <- colorRampPalette(c("#3C99B2","#ffffff","#EF2D00"), space = "rgb")(20)
pheatmap(t(metabolomics_noControls),
         color = scaleRYG,
         annotation_colors = my_colour,
         annotation_col = patient_annotation,
         annotation_row = metabolite_annotation,
         cluster_cols = T,
         #scale = "column",
         show_colnames = F,
         show_rownames = F,
         #scale = "row",
         main = "GC-MS Small Molecule COVID-19 HeatMap")


##### Exploratory Figures #####
metabolites_Known_Unknown <- t(metabolomics_noControls)
rownames(metabolites_Known_Unknown) <- metabolite_class$standardized_name
metabolites_Unknown <- metabolites_Known_Unknown[grep("unknown",rownames(metabolites_Known_Unknown)),]
metabolites_Known <- metabolites_Known_Unknown[-grep("unknown",rownames(metabolites_Known_Unknown)),]

# Boxplot of unknown metabolites
boxplot(t(metabolites_Unknown),
        main = "Boxplot of unknown GC-MS metabolites",
        ylab = "Log2(Intensity)")

par(mar = c(22.1, 4.1, 4.1, 4.1) # change the margins bottow, left, top, and right
    #lwd = 2, # increase the line thickness
    #cex.axis = 1.2 # increase default axis label size
)

# Boxplot of known metabolites
boxplot(t(metabolites_Known),
        main = "Boxplot of known GC-MS metabolites",
        ylab = "Log2(Intensity)",
        las = 1,
        xaxt = "n")

text(x = 1:length(rownames(metabolites_Known)),
     y = par("usr")[3] - 0.30, #subtract from y axis to push labels down
     labels = sub(" RT.*","",rownames(metabolites_Known)),
     xpd = NA, #print below axis
     srt = 90,
     adj = 1)


# Histogram of Unknowns and Knowns
par(mar = c(6.1, 4.1, 4.1, 4.1) # change the margins bottow, left, top, and right
    #lwd = 2, # increase the line thickness
    #cex.axis = 1.2 # increase default axis label size
)

hist(t(metabolites_Unknown),
     main = "Histogram of GC-MS Uknown features Log2(Intensity)",
     xlab = "Log2(Intensity)",
     xlim = c(0,30),
     las = 1,
     breaks = 20)

hist(t(metabolites_Known),
     main = "Histogram of GC-MS Known metabolites Log2(Intensity)",
     xlab = "Log2(Intensity)",
     xlim = c(0,35),
     las = 1,
     breaks = 20)

##### Fold Change to Median NONCOVID #####

metabolomics_COVID <- metabolomics_noControls[-grep("NONCOVID",rownames(metabolomics_noControls)),]
#write.csv(metabolomics_COVID, file = "H:/Projects/COVID19/SmallMolecule/Files/COVID_GCMetabolomics_normalized.csv")
metabolomics_NONCOVID <- metabolomics_noControls[grep("NONCOVID",rownames(metabolomics_noControls)),]
#write.csv(metabolomics_NONCOVID, file = "H:/Projects/COVID19/SmallMolecule/Files/NONCOVID_GCMetabolomics_normalized.csv")


par(mar = c(6.1, 4.1, 4.1, 4.1))
boxplot(t(metabolomics_NONCOVID), 
     main = "Boxplot of dynamic range for each Non-COVID patient",
     ylab = "Log2(Intensity",
     las=1,
     xaxt= "n")
text(x = 1:length(rownames(metabolomics_NONCOVID)),
     y = par("usr")[3] - 0.45,
     labels = rownames(metabolomics_NONCOVID),
     xpd = NA,
     srt = 35,
     adj = 1)

# Calculate the median for each column = metabolite/feature
NONCOVID_Median.MetaboliteIntensity <- colMedians(as.matrix(metabolomics_NONCOVID))

# Foldchange function -> take every row and divide it by the median vector
fold_change.FUNC <- function(x) x-NONCOVID_Median.MetaboliteIntensity
df_fc <- apply(metabolomics_COVID,1, fold_change.FUNC)

#n < -40
#palette <- distinctColorPalette(n)



my_colour = list(
  Gender = c(`F` = "#C0B9DD", `M` = "#234947"),
  ICU_1 = c(`No` = "#D0D3CA", `Yes` = "#368BA4"),
  Mech_Ventilation = c(`No` = "#D0D3CA", `Yes` = "#463F3A"))

patient_annotation$ICU_1 <- as.character(patient_annotation$ICU_1)
patient_annotation$ICU_1 <- sub("0","No",patient_annotation$ICU_1)
patient_annotation$ICU_1 <- sub("1","Yes",patient_annotation$ICU_1)
patient_annotation$ICU_1 <- factor(patient_annotation$ICU_1, levels = c("No","Yes"))

patient_annotation$Mech_Ventilation <- as.character(patient_annotation$Mech_Ventilation)
patient_annotation$Mech_Ventilation <- sub("0","No",patient_annotation$Mech_Ventilation)
patient_annotation$Mech_Ventilation <- sub("1","Yes",patient_annotation$Mech_Ventilation)
patient_annotation$Mech_Ventilation <- factor(patient_annotation$Mech_Ventilation, levels = c("No","Yes"))

# Heatmap of the Fold Change calculated from the median in NONCOVID cohort
par(mar = c(20, 4.1, 4.1, 4.1))
scaleRYG <- colorRampPalette(c("#3C99B2","#ffffff","#EF2D00"), space = "rgb")(20)


pheatmap(df_fc,
         color = scaleRYG,
         annotation_colors = my_colour,
         annotation_col = patient_annotation[-grep("NONCOVID",rownames(patient_annotation)),-c(1)],
         #annotation_row = metabolite_annotation,
         cluster_cols = T,
         #scale = "column",
         show_colnames = F,
         show_rownames = F,
         #scale = "row",
         main = "Heatmap Fold Change COVID/NONCOVID(median)")
  
par(mar = c(6.1, 4.1, 4.1, 4.1))
hist(df_fc,
     main = "Histogram of Fold Changes COVID relative to NONCOVID",
     xlab = "Fold Change of Log2",
     xlim = c(-15,15),
     breaks = 20)

# subset features with fold change greater than 1.2 to cluster
important_features_foldchange <- df_fc[rowSums(abs(df_fc)>2.2)>0,]

pheatmap(important_features_foldchange,
         color = scaleRYG,
         annotation_colors = my_colour,
         annotation_col = patient_annotation[-grep("NONCOVID",rownames(patient_annotation)),-c(1)],
         annotation_row = metabolite_annotation,
         cluster_cols = T,
         #scale = "column",
         show_colnames = F,
         show_rownames = F,
         #scale = "row",
         main = "Heatmap Fold Change COVID/NONCOVID(median)")

# Protein groups with fold change greater than 1.2
table(rowSums(abs(df_fc)>1.2)>0)

# Protein groups with fold change greater than 2
table(rowSums(abs(df_fc)>2)>0)

# Protein groups with fold change greater than 4
table(rowSums(abs(df_fc)>4)>0)

# Protein groups with fold change greater than 8
table(rowSums(abs(df_fc)>8)>0)


# bin protein groups so that we describe the changes within fold changes of -2 and 2
important_features_foldchange_compress <- important_features_foldchange
important_features_foldchange_compress[important_features_foldchange_compress < -2] = -2.1
important_features_foldchange_compress[important_features_foldchange_compress > 2] = 2.1

pheatmap(important_features_foldchange_compress,
         color = scaleRYG,
         annotation_colors = my_colour,
         annotation_col = patient_annotation[-grep("NONCOVID",rownames(patient_annotation)),-c(1)],
         annotation_row = metabolite_annotation,
         cluster_cols = T,
         #scale = "column",
         show_colnames = F,
         show_rownames = F,
         #scale = "row",
         main = "Heatmap Fold Change COVID/NONCOVID(median)")
