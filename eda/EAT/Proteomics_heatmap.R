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
out <- pheatmap(t(proteomics_noControls),
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


# Plot extracted gene-t0-cluster assignments
plot(out$tree_col)



#### Coagulation Cascade ####

coagulation_proteins_meta <- proteins_meta[which(proteins_meta$majority.protein %in% Deane_proteins$UniProt.ID),]

coagulation_proteomics <- proteomics_noControls[,which(colnames(proteomics_noControls) %in% coagulation_proteins_meta$biomolecule_id)]
colnames(coagulation_proteomics) <- coagulation_proteins_meta$majority.protein

coagulation_proteomics_COVID <- coagulation_proteomics[-grep("NONCOVID",rownames(coagulation_proteomics)),]
coagulation_proteomics_NONCOVID <- coagulation_proteomics[grep("NONCOVID",rownames(coagulation_proteomics)),]

out <- pheatmap(t(coagulation_proteomics),
         color = scaleRYG,
         annotation_colors = my_colour,
         annotation_col = patient_annotation,
         #annotation_row = metabolite_annotation,
         cluster_cols = T,
         #scale = "column",
         show_colnames = F,
         show_rownames = F,
         #scale = "row",
         main = "Coagulation Proteins of Interest\n COVID-19 patients only HeatMap")

plot(out$tree_row,
     main = "Dendogram of Proteins in Subset - coagulation ")
plot(out$tree_col,
     main = "Dendogram of Patients in subset - coagulation")

##### Exploratory Figures #####

proteomics_COVID <- t(proteomics_noControls)
rownames(proteomics_COVID) <- proteins_meta$majority.protein
proteomics_NONCOVID <- proteomics_COVID[,grep("NONCOVID",colnames(proteomics_COVID))]
proteomics_COVID <- proteomics_COVID[,-grep("NONCOVID",colnames(proteomics_COVID))]

# Boxplot of unknown proteomics
boxplot(proteomics_COVID,
        main = "Boxplot of COVID proteomics",
        ylab = "Log2(Intensity)")

par(mar = c(22.1, 4.1, 4.1, 4.1) # change the margins bottow, left, top, and right
    #lwd = 2, # increase the line thickness
    #cex.axis = 1.2 # increase default axis label size
)

# Boxplot of COVID proteomics
boxplot(proteomics_COVID,
        main = "Boxplot of COVID proteomics",
        ylab = "Log2(Intensity)",
        las = 1,
        xaxt = "n")

text(x = 1:length(colnames(proteomics_COVID)),
     y = par("usr")[3] - 0.30, #subtract from y axis to push labels down
     labels = colnames(proteomics_COVID),
     xpd = NA, #print below axis
     srt = 90,
     adj = 1)

# > median(proteomics_COVID)
# [1] 26.30916
# > mean(proteomics_COVID)
# [1] 26.72138

par(mar = c(10, 4.1, 4.1, 4.1) # change the margins bottow, left, top, and right
    #lwd = 2, # increase the line thickness
    #cex.axis = 1.2 # increase default axis label size
)
# Boxplot of NONCOVID proteomics
boxplot(proteomics_NONCOVID,
        main = "Boxplot of NONCOVID proteomics",
        ylab = "Log2(Intensity)",
        las = 1,
        xaxt = "n")

text(x = 1:length(colnames(proteomics_NONCOVID)),
     y = par("usr")[3] - 0.30, #subtract from y axis to push labels down
     labels = colnames(proteomics_NONCOVID),
     xpd = NA, #print below axis
     srt = 90,
     adj = 1)


# Histogram of Unknowns and Knowns
par(mar = c(6.1, 4.1, 4.1, 4.1) # change the margins bottow, left, top, and right
    #lwd = 2, # increase the line thickness
    #cex.axis = 1.2 # increase default axis label size
)

ICU <- patient_annotation[which(patient_annotation$ICU_1 == "Yes"),]
noICU <- patient_annotation[which(patient_annotation$ICU_1 == "No"),]

# NONCOVID ICU vs. noICU histograms
noICU_NONCOVID <- noICU[grep("NONCOVID",rownames(noICU)),]
ICU_NONCOVID <- ICU[grep("NONCOVID",rownames(ICU)),]
proteomics_NONCOVID_noICU <- proteomics_NONCOVID[,which(colnames(proteomics_NONCOVID) %in% rownames(noICU_NONCOVID))]
proteomics_NONCOVID_ICU <- proteomics_NONCOVID[,which(colnames(proteomics_NONCOVID) %in% rownames(ICU_NONCOVID))]

hist(proteomics_NONCOVID_ICU,
     main = "Histogram of NONCOVID in ICU samples\n features Log2(Intensity)",
     xlab = "Log2(Intensity)",
     #xlim = c(0,30),
     las = 1,
     breaks = 20)
abline(v= median(proteomics_NONCOVID_ICU), col = "red", lwd = 2)
# median of NONCOVID patients in the ICU 26.00141

hist(proteomics_NONCOVID_noICU,
     main = "Histogram of NONCOVID NOT in the ICU samples\n features Log2(Intensity)",
     xlab = "Log2(Intensity)",
     #xlim = c(0,30),
     las = 1,
     breaks = 20)
abline(v= median(proteomics_NONCOVID_noICU), col = "red", lwd = 2)
# median of NONCOVID patients in the ICU 25.95939

# COVID ICU vs noICU histograms
noICU_COVID <- noICU[-grep("NONCOVID",rownames(noICU)),]
ICU_COVID <- ICU[-grep("NONCOVID",rownames(ICU)),]
proteomics_COVID_noICU <- proteomics_COVID[,which(colnames(proteomics_COVID) %in% rownames(noICU_COVID))]
proteomics_COVID_ICU <- proteomics_COVID[,which(colnames(proteomics_COVID) %in% rownames(ICU_COVID))]

hist(proteomics_COVID_ICU,
     main = "Histogram of COVID in ICU samples\n features Log2(Intensity)",
     xlab = "Log2(Intensity)",
     #xlim = c(0,30),
     las = 1,
     breaks = 20)
abline(v= median(proteomics_COVID_ICU), col = "red", lwd = 2)
# median of NONCOVID patients in the ICU 26.38699

hist(proteomics_COVID_noICU,
     main = "Histogram of COVID NOT in the ICU samples\n features Log2(Intensity)",
     xlab = "Log2(Intensity)",
     #xlim = c(0,30),
     las = 1,
     breaks = 20)
abline(v= median(proteomics_COVID_noICU), col = "red", lwd = 2)
# median of NONCOVID patients in the ICU 26.22727

##### Fold Change to Median NONCOVID #####

# Calculate the median for each column = metabolite/feature
NONCOVID_Median.ProteinIntensity <- rowMedians(as.matrix(proteomics_NONCOVID))

# Foldchange function -> take every row and divide it by the median vector
fold_change.FUNC <- function(x) x-NONCOVID_Median.ProteinIntensity
df_fc <- apply(proteomics_COVID,2, fold_change.FUNC)

# New color palette patients
my_colour = list(
  Gender = c(`F` = "#C0B9DD", `M` = "#234947"),
  Batch = c(`1` = "#FDF0FE", `2` = "#E8C5E9", `3` = "#9FA0C3", `4`  = "#8B687F", `5` = "#7B435B", `6` = "#5C3344", `7` = "#1C093D"),
  ICU_1 = c(`No` = "#D0D3CA", `Yes` = "#368BA4"),
  Mech_Ventilation = c(`No` = "#D0D3CA", `Yes` = "#463F3A"))

# Heatmap of the Fold Change calculated from the median in NONCOVID cohort
scaleRYG <- colorRampPalette(c("#3C99B2","#ffffff","#EF2D00"), space = "rgb")(20)

out <- pheatmap(df_fc,
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

plot(out$tree_col,
     main = "Fold Change COVID/NONCOVID(median)\n Patients")


#Export fold change data frame ordered by protein clustering
# Re-0rder original data (proteins) to match ordering in heatmap (top-to-bottom)
df_fc_export <- df_fc[rownames(df_fc[out$tree_row[["order"]],]),colnames(df_fc[,out$tree_col[["order"]]])]
write.csv(df_fc_export, file = "P:/All_20200428_COVID_plasma_multiomics/Proteomics/EAT_unsupervised_analysis/COVID_proteomics_Heatmap_COVIDrelativeNONCOVIDmedian.csv")

# Histogram of Fold Changes
par(mar = c(6.1, 4.1, 4.1, 4.1))
hist(df_fc,
     main = "Histogram of Fold Changes COVID relative to NONCOVID",
     xlab = "Fold Change in Log2 space",
     xlim = c(-15,15),
     breaks = 50)

# Protein groups with fold change greater than 1.2
table(rowSums(abs(df_fc)>1.2)>0)

# Protein groups with fold change greater than 2
table(rowSums(abs(df_fc)>2)>0)

# Protein groups with fold change greater than 4
table(rowSums(abs(df_fc)>4)>0)

# Protein groups with fold change greater than 8
table(rowSums(abs(df_fc)>8)>0)


# subset features with fold change greater than 4.2 to cluster
important_features_foldchange <- df_fc[rowSums(abs(df_fc)>6.2)>0,]

out <- pheatmap(important_features_foldchange,
                color = scaleRYG,
                annotation_colors = my_colour,
                annotation_col = patient_annotation[-grep("NONCOVID",rownames(patient_annotation)),-c(1)],
                #annotation_row = metabolite_annotation,
                cluster_cols = T,
                #scale = "column",
                show_colnames = F,
                show_rownames = F,
                #scale = "row",
                main = "Heatmap Fold Change COVID/NONCOVID(median)\n Subset Features that have abs(FC) > 6.2 in at least one patient")

##### Fold Change Median of NONCOVID not in ICU ####
# Calculate the median for each column = metabolite/feature
NONCOVID_NOICU_Median.ProteinIntensity <- rowMedians(as.matrix(proteomics_NONCOVID_noICU))

# Foldchange function -> take every row and divide it by the median vector
fold_change.FUNC <- function(x) x-NONCOVID_NOICU_Median.ProteinIntensity
df_fc <- apply(proteomics_COVID,2, fold_change.FUNC)

# Heatmap of the Fold Change calculated from the median in NONCOVID cohort
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
         main = "Heatmap Fold Change COVID/NONCOVID not in the ICU(median)")

##### Fold Change Median of NONCOVID not in ICU ####
# Calculate the median for each column = metabolite/feature
NONCOVID_NOICU_Mean <- mean(proteomics_NONCOVID_noICU)

# Foldchange function -> take every row and divide it by the median vector
fold_change.FUNC <- function(x) x-NONCOVID_NOICU_Mean
df_fc <- apply(coagulation_proteomics_COVID,1, fold_change.FUNC)

# Heatmap of the Fold Change calculated from the median in NONCOVID cohort
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
         main = "Heatmap Fold Change COVID/NONCOVID not in the ICU(mean)")

##### Fold Change in Coagulation Subset ######

# Calculate the median for each column = metabolite/feature
NONCOVID_Coagulation_Median.ProteinIntensity <- colMedians(as.matrix(coagulation_proteomics_NONCOVID))

# Foldchange function -> take every row and divide it by the median vector
fold_change.FUNC <- function(x) x-NONCOVID_Coagulation_Median.ProteinIntensity
df_fc <- apply(coagulation_proteomics_COVID,1, fold_change.FUNC)

# Heatmap of the Fold Change calculated from the median in NONCOVID cohort
scaleRYG <- colorRampPalette(c("#3C99B2","#ffffff","#EF2D00"), space = "rgb")(20)

out <- pheatmap(df_fc,
         color = scaleRYG,
         annotation_colors = my_colour,
         annotation_col = patient_annotation[-grep("NONCOVID",rownames(patient_annotation)),-c(1)],
         #annotation_row = metabolite_annotation,
         cluster_cols = T,
         #scale = "column",
         show_colnames = F,
         show_rownames = F,
         #scale = "row",
         main = "Heatmap Fold Change COVID/NONCOVID(median) in Coagulation Subset Proteins")

plot(out$tree_row,
     main = "Dendogram of Proteins\n Fold Change COVID/NONCOVID(median) in Coagulation Subset")

plot(out$tree_col,
     main = "Dendogram of COVID Patients\n Fold Change COVID/NONCOVID(median) in Coagulation Subset")
