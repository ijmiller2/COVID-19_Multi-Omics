##### 19_KAO_Clinical_metrics_correlations.R #######

library(DBI)
library(RSQLite)

source("eda/KAO/0_pathway_toolkit.R")
source("eda/KAO/0_two_omic_cor_v2.R")
library(pheatmap)
library(igraph)


colors <- read.csv("reference/color_palette.txt",  stringsAsFactors = F)[,2]
palette(colors)


### Fig 2A - PCA #### 

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

clinical_data <- dbGetQuery(con, "SELECT *
                            FROM deidentified_patient_metadata")

df_metabolites<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, metabolomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45, WHO_ordinal_at_day_28
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")

df_lipids<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, lipidomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45, WHO_ordinal_at_day_28
           FROM lipidomics_measurements
           INNER JOIN lipidomics_runs ON lipidomics_runs.replicate_id = lipidomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = lipidomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = lipidomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_proteins<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, proteomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45, WHO_ordinal_at_day_28
           FROM proteomics_measurements
           INNER JOIN proteomics_runs ON proteomics_runs.replicate_id = proteomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = proteomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = proteomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_transcripts<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, transcriptomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45, WHO_ordinal_at_day_28
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

#### Creating a wide-format data frame  #####

df_wide <- reshape(df, timevar = "biomolecule_id", v.names = "normalized_abundance",
                   idvar = "sample_id", direction = "wide" )

names(df_wide)
df_wide_all <- df_wide[!(df_wide$sample_id == 88 | df_wide$sample_id == 81 | df_wide$sample_id == 129), ]

df_wide_all$sample_id
clinical_data$sample_id

merge_df <- merge(clinical_data, df_wide_all, by ="sample_id")

#### correlation between clinical data and biomolecules #### 
clinical_exprs <- t(merge_df[,c(13:18, 24:32)])
df_exprs <- t(merge_df[, -c(1:36)])
dim(clinical_exprs)
dim(df_exprs)

clinic_cor <- two_omic_cor_fast(df_exprs, clinical_exprs, method = "kendall")

#### biomolecules_information from db #### 

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

df_proteins<- dbGetQuery(con, "SELECT * FROM metadata
            WHERE metadata_type = 'fasta_header'  
           ")

biomolecules <-dbGetQuery(con, "SELECT * FROM biomolecules")

pvalues <- dbGetQuery(con, "SELECT * FROM pvalues WHERE comparison = 'COVID_vs_NONCOVID'")
pvalues_HFD <- dbGetQuery(con, "SELECT * FROM pvalues WHERE comparison = 'Hospital_free_days_45'")


# pull
metadata_GO_bp <- dbGetQuery(con, "SELECT * FROM metadata WHERE metadata_type = 'GO_biological_process'")
metadata_GO_mf <- dbGetQuery(con, "SELECT * FROM metadata WHERE metadata_type = 'GO_molecular_function'")
metadata_class <- dbGetQuery(con, "SELECT * FROM metadata WHERE metadata_type = 'Larger_class' OR metadata_type = 'Lipid Class'")

dbDisconnect(con)

#### biomolecules annotation #### 

df_proteins$geneNames<- apply(df_proteins , 1, function(x) strsplit(strsplit(x[4], "GN=")[[1]][2], " ")[[1]][1])
biomolecules_merge <- merge(biomolecules, df_proteins[,c(2,4,5)], by = "biomolecule_id", all.x = T)
biomolecules_df <- merge(biomolecules_merge, metadata_GO_bp[c(2,4)], by = "biomolecule_id", all.x = T)

names(biomolecules_df)[5] <- "fasta_header"
names(biomolecules_df)[7] <- "GO_BP"

biomolecules_df$standardized_name[biomolecules_df$omics_id == 1] <- biomolecules_df$geneNames[biomolecules_df$omics_id == 1]

#### which clinical - feature correlations are > 0.4 #### 

hist(clinic_cor)
filter_clinical <- rowSums(clinic_cor > 0.3) > 0
table(filter_clinical)


filter_features <- colSums(clinic_cor > 0.3) > 0
table(filter_features)

dim(clinic_cor)
pheatmap(clinic_cor[ ,filter_features])

### write out RData file #### 

colnames(clinic_cor) <-  sub("normalized_abundance.", "",colnames(clinic_cor))

save(clinic_cor, biomolecules_df, file = "data/Correlation_with_clinical_measurements_KAO.RData")
clinic_cor_std <- clinic_cor
colnames(clinic_cor_std) <- biomolecules_df$standardized_name[match(colnames(clinic_cor), biomolecules_df$biomolecule_id)]

write.csv(clinic_cor_std[,filter_features], "clinical_correlation.csv")

### GO terms for WBC counts ####

class_and_GO_bp <- make_reference_sets(biomolecules_df$GO_BP, as.character(biomolecules_df$biomolecule_id))

row.names(clinic_cor)
biomolecules_df$standardized_name[ as.character(biomolecules_df$biomolecule_id) %in% colnames(clinic_cor)]


#enrich D-Dimer
enrichment_ddimer <- enrichment(colnames(clinic_cor)[clinic_cor[3, ] > 0.3], class_and_GO_bp, colnames(clinic_cor))

## procalcitonin 

enrichment_procalcitonin <- enrichment(colnames(clinic_cor)[clinic_cor[4, ] > 0.3], class_and_GO_bp, colnames(clinic_cor))

## white blood cell count 
enrichment_wbc <- enrichment(colnames(clinic_cor)[clinic_cor[8, ] > 0.3], class_and_GO_bp, colnames(clinic_cor))

pdf("plots/Enrichment_WBC_percent_correlation.pdf", width = 7, height = 2.5)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)

barplot(-log10(enrichment_wbc [order(enrichment_wbc[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 1, 
        names = enrichment_wbc [order(enrichment_wbc[,2]),1][1:5][5:1],
        main = "GO Biological processes WBC")
dev.off()


## Neutrophils 
enrichment_neutrophil <- enrichment(colnames(clinic_cor)[clinic_cor[12, ] > 0.30], class_and_GO_bp, colnames(clinic_cor))

pdf("plots/Enrichment_neutrophil_percent_correlation.pdf", width = 7, height = 2.5)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)

barplot(-log10(enrichment_neutrophil [order(enrichment_neutrophil [,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 1, 
        bty = "n",
        names = enrichment_neutrophil [order(enrichment_neutrophil[,2]),1][1:5][5:1],
        main = "GO Biological processes Neutrophil")
dev.off()

## Lymphocytes 
enrichment_lymphocytes <- enrichment(colnames(clinic_cor)[clinic_cor[13, ] > 0.3], class_and_GO_bp, colnames(clinic_cor))

pdf("plots/Enrichment_lymphocyte_percent_correlation.pdf", width = 7, height = 2.5)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)

barplot(-log10(enrichment_lymphocytes [order(enrichment_lymphocytes [,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 1, 
        names = enrichment_lymphocytes [order(enrichment_lymphocytes[,2]),1][1:5][5:1],
        main = "GO Biological processes Lymphocytes")
dev.off()


## Monocyte 
enrichment_monocytes <- enrichment(colnames(clinic_cor)[clinic_cor[14, ] > 0.30], class_and_GO_bp, colnames(clinic_cor))

pdf("plots/Enrichment_monocyte_percent_correlation.pdf", width = 7, height = 2.5)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)

barplot(-log10(enrichment_monocytes [order(enrichment_monocytes [,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 1, 
        names = enrichment_monocytes [order(enrichment_monocytes[,2]),1][1:5][5:1],
        main = "GO Biological processes Monocytes")
dev.off()


## Eosinophils
enrichment_eosinophils <- enrichment(colnames(clinic_cor)[clinic_cor[15, ] > 0.30], class_and_GO_bp, colnames(clinic_cor))

pdf("plots/Enrichment_eosinophil_percent_correlation.pdf", width = 7, height = 2.5)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)

barplot(-log10(enrichment_eosinophils [order(enrichment_eosinophils [,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 1, 
        names = enrichment_eosinophils [order(enrichment_eosinophils[,2]),1][1:5][5:1],
        main = "GO Biological processes Eosinophils")
dev.off()


## platelets 
enrichment_platelets <- enrichment(colnames(clinic_cor)[clinic_cor[11, ] > 0.30], class_and_GO_bp, colnames(clinic_cor))

pdf("plots/Enrichment_platelet_percent_correlation.pdf", width = 7, height = 2.5)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)

barplot(-log10(enrichment_platelets [order(enrichment_platelets[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 1, 
        names = enrichment_platelets [order(enrichment_platelets[,2]),1][1:5][5:1],
        main = "GO Biological processes Platelets")
dev.off()

