### 18_KAO_Comparing_WHO_score_to_HFD.R ###


## These data pull infromation from pvalues table and provides numbers for the text about figure 2. 

library(DBI)
library(RSQLite)

colors <- read.csv("reference/color_palette.txt",  stringsAsFactors = F)[,2]
palette(colors)

##### Connecting to the DB #### 

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")


df_metabolites <- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, metabolomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45, WHO_ordinal_at_day_28, SOFA, SAPSII
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")

df_lipids<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, lipidomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45, WHO_ordinal_at_day_28, SOFA, SAPSII
           FROM lipidomics_measurements
           INNER JOIN lipidomics_runs ON lipidomics_runs.replicate_id = lipidomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = lipidomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = lipidomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_proteins<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, proteomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45, WHO_ordinal_at_day_28, SOFA, SAPSII
           FROM proteomics_measurements
           INNER JOIN proteomics_runs ON proteomics_runs.replicate_id = proteomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = proteomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = proteomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_transcripts<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, transcriptomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45, WHO_ordinal_at_day_28, SOFA, SAPSII
           FROM transcriptomics_measurements
           INNER JOIN transcriptomics_runs ON transcriptomics_runs.replicate_id = transcriptomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = transcriptomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = transcriptomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")
dbDisconnect(con)

df <- rbind(df_metabolites, df_lipids, df_proteins, df_transcripts)

df <- df[df$sample_id != 54, ]
names(df)

#### Creating a wide-format data frame #####

df_wide <- reshape(df, timevar = "biomolecule_id", v.names = "normalized_abundance",
                   idvar = "sample_id", direction = "wide" )

names(df_wide)

##### Plotting HFD vs. WHO score ##### 

plot(df_wide$Hospital_free_days_45, df_wide$WHO_ordinal_at_day_28)

df_COVID <- df_wide[ df_wide$COVID == 1, 1:7]
df_nonCOVID <- df_wide[df_wide$COVID == 0, 1:7]

pdf("E://COVID-19_Multi-Omics/plots/HFD45_vs_WHOscore_v1.pdf", width = 5, height = 5, useDingbats = F)
plot(df_COVID$Hospital_free_days_45[order(df_COVID$WHO_ordinal_at_day_28)], pch = 19, col = 3, ylab = "Hospital Free Days 45" , las = 1,
     main = "Hospital free days vs. WHO ordinal score")
par(new = T)
plot(df_COVID$WHO_ordinal_at_day_28[order(df_COVID$WHO_ordinal_at_day_28)], pch = 19, col = 6, ylim = c(0,8), axes =  F, ylab = "")
axis(4, at = c(0,2,4,6,8), las = 1)
dev.off()


pdf("E://COVID-19_Multi-Omics/plots/HFD45_vs_WHOscore.pdf", width = 5, height = 5, useDingbats = F)
plot(df_COVID$Hospital_free_days_45 ~ jitter(df_COVID$WHO_ordinal_at_day_28-0.3, 0.4), xlim=c(-0.5,8.5),
     xlab = "WHO ordinal at day 28", col = 3, ylab = "Hospital free days 45", main = "Comparing Hospital free days\n and WHO ordinal score")
par(new = T)
x <- boxplot(df_COVID$Hospital_free_days_45 ~ df_COVID$WHO_ordinal_at_day_28, boxwex = 0.2, col = 3,
             at = c(0.5,3.5,4.5,5.5,6.5,8.5), new = T, xaxt = "n", yaxt  = "n", outline = F)

dev.off()

hist(df_COVID$Hospital_free_days_45, breaks = 20, col = 3,
     main = "Histogram\n Hospital Free Days-45",
     xlab = "Hospital free days-45")
hist(df_COVID$WHO_ordinal_at_day_28, breaks = 20, col = 6, 
     main = "Histogram\n WHO ordinal at day 28", 
     xlab = "WHO ordinal at day 28")

plot(df_COVID$Hospital_free_days_45, df_COVID$WHO_ordinal_at_day_28)
plot(df_nonCOVID$Hospital_free_days_45, df_nonCOVID$WHO_ordinal_at_day_28)
hist(df_nonCOVID$Hospital_free_days_45, breaks = 20)

## Creating color scale for HFD
breaks = seq(0,45, by =5)
bin_values <- .bincode(df_wide_all$Hospital_free_days_45, breaks, include.lowest = T)
col_HFD<- colorRampPalette(c("black", "gray75"))(9)[bin_values]
