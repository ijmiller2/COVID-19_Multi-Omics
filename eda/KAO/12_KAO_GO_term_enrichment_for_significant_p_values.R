#### 12_KAO_GO_term_enrichment_for_significant_p_values.R #### 

## testing GO term enrichment of biomolecules significant for COVID 

library(DBI)
library(RSQLite)

source("eda/KAO/0_pathway_toolkit.R")
library(pheatmap)
colors <- read.csv("reference/color_palette.txt",  stringsAsFactors = F)[,2]
palette(colors)

##### getting go terms and pvalues and biomoleucle class ##### 

## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# pull
metadata_GO_bp <- dbGetQuery(con, "SELECT * FROM metadata WHERE metadata_type = 'GO_biological_process'")
metadata_GO_mf <- dbGetQuery(con, "SELECT * FROM metadata WHERE metadata_type = 'GO_molecular_function'")
metadata_class <- dbGetQuery(con, "SELECT * FROM metadata WHERE metadata_type = 'Larger_class' OR metadata_type = 'Lipid Class'")

pvalues <- dbGetQuery(con, "SELECT * FROM pvalues WHERE comparison = 'COVID_vs_NONCOVID'")
pvalues_HFD <- dbGetQuery(con, "SELECT * FROM pvalues WHERE comparison = 'Hospital_free_days_45'")

df_metabolites<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, metabolomics_measurements.biomolecule_id, COVID, Age_less_than_90, Gender, ICU_1, Hospital_free_days_45, Charlson_score, SOFA
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")



df_lipids<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, lipidomics_measurements.biomolecule_id, COVID, Age_less_than_90, Gender, ICU_1, Hospital_free_days_45, Charlson_score, SOFA
           FROM lipidomics_measurements
           INNER JOIN lipidomics_runs ON lipidomics_runs.replicate_id = lipidomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = lipidomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = lipidomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_proteins<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, proteomics_measurements.biomolecule_id, COVID, Age_less_than_90, Gender, ICU_1, Hospital_free_days_45, Charlson_score, SOFA
           FROM proteomics_measurements
           INNER JOIN proteomics_runs ON proteomics_runs.replicate_id = proteomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = proteomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = proteomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")

# disconnect
dbDisconnect(con) 

df <- rbind(df_metabolites, df_lipids, df_proteins)

df <- df[df$sample_id != 54, ]

##### make reference set with go term bp and class info #### 
metadata_GO_bp_class<- rbind(metadata_GO_bp, metadata_class)
metadata_GO_mf_class <- rbind(metadata_GO_mf, metadata_class)

class_and_GO_bp <- make_reference_sets(metadata_GO_bp_class$metadata_value, metadata_GO_bp_class$biomolecule_id)
class_and_GO_mf <- make_reference_sets(metadata_GO_mf_class$metadata_value, metadata_GO_mf_class$biomolecule_id)


#### performing enrichment based LR test significant for COVID #### 
enrichment_significant_pvalues <- enrichment(pvalues$biomolecule_id[pvalues$q_value < 0.05], class_and_GO_bp, pvalues$biomolecule_id)
enrichment_mf_sig_pvalues <- enrichment(pvalues$biomolecule_id[pvalues$q_value < 0.05], class_and_GO_mf, pvalues$biomolecule_id)

pdf("12_KAO_GO_BP_Enrichment_based_on_LRtest_COVID.pdf")
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
barplot(-log10(enrichment_significant_pvalues[order(enrichment_significant_pvalues[,2]),2][1:10][10:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 5, 
        names = enrichment_significant_pvalues[order(enrichment_significant_pvalues[,2]),1][1:10][10:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significant with COVID")
dev.off()

#barplot(-log10(enrichment_mf_sig_pvalues[order(enrichment_mf_sig_pvalues[,2]),2][1:10][10:1]), xlab ="-log10(p-value)", cex.names = 0.9, horiz = T, col = 5, names = enrichment_mf_sig_pvalues[order(enrichment_mf_sig_pvalues[,2]),1][1:10][10:1], main = "GO_MF enrichment\nbiomolecules significant with COVID")

#### Extracting out change in mean based on COVID status using summary(lm())$coefficients
lm_coefficients <- function(biomolecule_id, formula_test, data, return_first = T, return_category = NULL){
  lm_formula_test <- lm(formula_test, data = data[data$biomolecule_id == biomolecule_id, ])
  if(return_first == T){
    c(coeficient =  summary(lm_formula_test)$coefficients[2,3], pvalue = summary(lm_formula_test)$coefficients[2,4])
  } else if (!is.null(return_category) & return_category %in% row.names(summary(lm_formula_test)$coefficients)){
    index <- match(return_category, row.names(summary(lm_formula_test)$coefficients))
    c(coeficient =  summary(lm_formula_test)$coefficients[index,3], pvalue = summary(lm_formula_test)$coefficients[index,4])
  }
}


coefficients_covid <- apply(pvalues, 1, function(x) 
  lm_coefficients(as.numeric(x[2]), formula_test = normalized_abundance ~ COVID + ICU_1 + Gender + Age_less_than_90, 
                  data = df))

up_covid <- coefficients_covid[1,] >0.5
down_covid <- coefficients_covid[1,] < -0.5

#### performing enrichment based LR test significant UP for COVID #### 
enrichment_up_significant_pvalues <- enrichment(pvalues$biomolecule_id[pvalues$q_value < 0.05 & up_covid], class_and_GO_bp, pvalues$biomolecule_id)
enrichment_down_significant_pvalues <-  enrichment(pvalues$biomolecule_id[pvalues$q_value < 0.05 & down_covid], class_and_GO_bp, pvalues$biomolecule_id)

pdf("plots/12_KAO_GO_BP_Enrichment_based_on_LRtest_UP_COVID.pdf")
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
par(mfrow = c(2,1))
barplot(-log10(enrichment_up_significant_pvalues[order(enrichment_up_significant_pvalues[,2]),2][1:10][10:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 5, 
        names = enrichment_up_significant_pvalues[order(enrichment_up_significant_pvalues[,2]),1][1:10][10:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significantly up with COVID")


barplot(-log10(enrichment_down_significant_pvalues[order(enrichment_down_significant_pvalues[,2]),2][1:10][10:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 6, 
        names = enrichment_down_significant_pvalues[order(enrichment_down_significant_pvalues[,2]),1][1:10][10:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significantly down with COVID")

dev.off()

#barplot(-log10(enrichment_mf_sig_pvalues[order(enrichment_mf_sig_pvalues[,2]),2][1:10][10:1]), xlab ="-log10(p-value)", cex.names = 0.9, horiz = T, col = 5, names = enrichment_mf_sig_pvalues[order(enrichment_mf_sig_pvalues[,2]),1][1:10][10:1], main = "GO_MF enrichment\nbiomolecules significant with COVID")

library(pheatmap)

df_wide <- reshape(df, timevar = "biomolecule_id", v.names = "normalized_abundance",
                   idvar = "sample_id", direction = "wide" )

names(df_wide)
df_wide_exprs <- df_wide[,-c(1:8)]

annotation_col <- as.data.frame(df_wide[,2:8])
row.names(annotation_col) = row.names(df_wide)


pheatmap(t(df_wide_exprs[,pvalues$q_value <0.05 & up_covid]), 
         scale = "row", annotation_col= annotation_col[c(1,4,5)], 
         breaks = c(-10, -4, -3, -2, -1, 0, 1, 2, 3, 4, 10),
         col = colorRampPalette(c(3,"white", 5))(10),
         show_rownames = F, show_colnames = F)


######## Enrichment HFD ####### 

enrichment_HFD <- enrichment(pvalues_HFD$biomolecule_id[pvalues_HFD$q_value < 0.05], class_and_GO_bp, pvalues_HFD$biomolecule_id)

pdf("12_KAO_GO_BP_Enrichment_based_on_lm_HFD.pdf")
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
barplot(-log10(enrichment_HFD[order(enrichment_HFD[,2]),2][1:10][10:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 5, 
        names = enrichment_HFD[order(enrichment_HFD[,2]),1][1:10][10:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significant with HFD")
dev.off()

######## Enrichment HFD _up_down #### 
coefficients_hfd <- apply(pvalues, 1, function(x) 
  lm_coefficients(as.numeric(x[2]), formula_test = normalized_abundance ~ Hospital_free_days_45, 
                  data = df))

hist(coefficients_hfd[1,])
up_hfd <- coefficients_hfd[1,] >  1
down_hfd <- coefficients_hfd[1,] < -1

enrichment_up_hfd <- enrichment(pvalues_HFD$biomolecule_id[pvalues_HFD$q_value < 0.05 & up_hfd], class_and_GO_bp, pvalues_HFD$biomolecule_id)
enrichment_down_hfd <-  enrichment(pvalues_HFD$biomolecule_id[pvalues_HFD$q_value < 0.05 & down_hfd], class_and_GO_bp, pvalues_HFD$biomolecule_id)

pdf("plots/12_KAO_GO_BP_Enrichment_based_on_Linear_regression_HFD.pdf")
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
par(mfrow = c(2,1))
barplot(-log10(enrichment_up_hfd[order(enrichment_up_hfd[,2]),2][1:10][10:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 5, 
        names = enrichment_up_hfd[order(enrichment_up_hfd[,2]),1][1:10][10:1],
        main = "Enrichment of molecule class and GO Biological processes\npositively correlated with HFD")


barplot(-log10(enrichment_down_hfd[order(enrichment_down_hfd[,2]),2][1:10][10:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 6, 
        names = enrichment_down_hfd[order(enrichment_down_hfd[,2]),1][1:10][10:1],
        main = "Enrichment of molecule class and GO Biological processes\nnegatively correlated with HFD")

dev.off()

###### Heatmap of sig with HFD ########

pheatmap(t(df_wide_exprs[,pvalues_HFD$q_value <0.05]), 
         scale = "row", annotation_col= annotation_col[c(1,4,5)], 
         breaks = c(-10, -4, -3, -2, -1, 0, 1, 2, 3, 4, 10),
         col = colorRampPalette(c(3,"white", 5))(10),
         show_rownames = F, show_colnames = F)

