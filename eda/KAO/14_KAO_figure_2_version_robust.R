##### 14_KAO_figure_2_version_robust.R ##### 

## These data pull infromation from pvalues table and provides numbers for the text about figure 2. 

library(DBI)
library(RSQLite)

source("eda/KAO/0_pathway_toolkit.R")
library(pheatmap)
library(ellipse)
library(scales)

colors <- read.csv("reference/color_palette.txt",  stringsAsFactors = F)[,2]
palette(colors)


### Fig 2A - PCA #### 

con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")


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

#### Creating a wide-format data frame to facilitate PCA #####

df_wide <- reshape(df, timevar = "biomolecule_id", v.names = "normalized_abundance",
                   idvar = "sample_id", direction = "wide" )

names(df_wide)
df_wide_all <- df_wide[!(df_wide$sample_id == 88 | df_wide$sample_id == 81 | df_wide$sample_id == 129), ]

df_wide_exprs <- df_wide_all[,-c(1:5)]


##### Priciple component analysis #######
pca1 <- prcomp(df_wide_exprs)


## Creating color scale for HFD
breaks = seq(0,45, by =5)
bin_values <- .bincode(df_wide_all$Hospital_free_days_45, breaks, include.lowest = T)
col_HFD<- colorRampPalette(c("black", "gray75"))(9)[bin_values]

## Creating color scale for WHO
bin_values_WHO <- df_wide_all$WHO_ordinal_at_day_28
col_WHO<- colorRampPalette(c("gray60", "black"))(9)[bin_values_WHO+1]

## creating ellipse with 95% CI for top and bottom bin values
lvs2 <- levels(as.factor(bin_values_WHO))
pts.array2 <- array(0, dim = c(100, 2, length(lvs2)))

for(i in c(1,7)) {
  inx <- bin_values_WHO == lvs2[i]
  groupVar <- var(cbind(pca1$x[inx, 1], pca1$x[inx,2]), na.rm = T )
  groupMean <- cbind(mean(pca1$x[inx, 1], na.rm= T), mean(pca1$x[inx, 2], na.rm =T))
  pts.array2[,,i] <- ellipse(groupVar, centre = groupMean, level = 0.95, npoints = 100)
}

## Plotting PCA 
pdf("plots/Fig2_PCA_all_omes_WHO_color.pdf", width = 5, height = 5, useDingbats = F)
plt_original <- par()$plt
par(mgp = c(1.8, 0.7, 0), tcl = -0.3)
par("plt" = c(0.92,0.97, 0.22, 0.43))
plot(breaks, type ="n", bty = "n", xaxt = "n", xlab = "", ylab = "WHO at day 28", las =1 , )
for(i in 1:length(breaks[-1])){
  rect(0,breaks[i],10,breaks[i+1], col = colorRampPalette(c("gray75", "black"))(9)[i], border =NA)
}

par(new = T)
par("plt" = plt_original)

plot(pca1$x, col= col_WHO, pch =c(17,19)[df_wide_all$COVID+1], cex = 1.5, lwd = 3, bty = "l",
     ylim = c(-125, 125), xlim = c(-175, 150), las =1 ,
     main = "Pricipal component anlaysis - all omes",
     xlab = paste("PC1 (", round(summary(pca1)$importance[2,1]*100, digits = 2), "%)", sep = ""),
     ylab = paste("PC2 (", round(summary(pca1)$importance[2,2]*100, digits = 2), "%)", sep = ""))


polygon(pts.array2[,,1], border = "gray75", lty = 2, lwd = 2)
polygon(pts.array2[,,7], border = "black", lty = 2, lwd = 2)

legend("bottomleft", pch = c(19, 17), c("COVID-19", "non-COVID-19"), bty = "n")


dev.off()

### Fig 2B - volcano plot COVID-status comparing robust vs. std. linear regression #### 

# connect 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# pull pvalues 

pvalues_COVID_robust <- dbGetQuery(con, "SELECT * FROM pvalues 
          INNER JOIN biomolecules on biomolecules.biomolecule_id = pvalues.biomolecule_id
          WHERE pvalues.comparison = 'COVID' AND pvalues.test = 'LR_test_robust' AND biomolecules.keep = 1
                  ")

# disconnect
dbDisconnect(con) 

## load in log2 fold change values 

FC <- read.csv("data/COVID_fc_by_biomolecule_ID.csv")

## load in transcript pvalues 

transcript_pvalues <- read.delim("P:/All_20200428_COVID_plasma_multiomics/Transcriptomics/AllGenes.txt", sep = " ")
transcript_pvalues <- data.frame(standardized_name = row.names(transcript_pvalues), PP = 1.0000000000001-transcript_pvalues[,1], stringsAsFactors = F)

## merging pvalues from database with DESeq results
pvalues_COVID_robust <- pvalues_COVID_robust[, -grep('biomolecule_id', names(pvalues_COVID_robust))[2]]

pvalues_COVID_robust_deseq <- merge(pvalues_COVID_robust, transcript_pvalues, by = "standardized_name", all.x = T)

## replacing original q-values based on log-likelyhood with DEseq results (1-PP)

pvalues_COVID_robust_deseq$q_value[!is.na(pvalues_COVID_robust_deseq$PP)] <- pvalues_COVID_robust_deseq$PP[!is.na(pvalues_COVID_robust_deseq$PP)] 


## merge pvalues and FC by biomolecule_id

m_robust <- merge(pvalues_COVID_robust_deseq, FC, by.x = 'biomolecule_id', by.y = 'biomolecule_id')
m_robust$omics_id <- as.factor(m_robust$omics_id)
m_robust$omics_id <- relevel(m_robust$omics_id, ref = 5)
m_robust <- m_robust[order(m_robust$omics_id), ]



unknown <- grepl("nknown", m_robust$standardized_name)
table(unknown)

## Plot Volcano plot 
pdf("plots/KAO_fig_2A_Volcano_plot_v2_robust_linear_regression.pdf", useDingbats = F, width = 5, height = 5)

par(mgp = c(1.8, 0.7, 0), tcl = -0.3)
plot(m_robust$FC[unknown], -log(m_robust$q_value)[unknown],
     pch = 1, col = "gray40",
     ylim = c(0, 34), xlim = c(-4, 6),
     main = "Effect of COVID", 
     ylab = "-log(FDR)",
     xlab = "log2(COVID-19/non-COVID-19)",
     las = 1, bty = "l", cex.axis = 1.2)
abline(h = -log(0.05), lty = 2)

points(m_robust$FC[!unknown], -log(m_robust$q_value)[!unknown], 
       pch = 19, col = colors[c(12,9,10,11,11)][m_robust$omics_id[!unknown]])

legend("topleft", c("transcripts", "proteins", "lipids", "small molecules", "unidentified features"), 
       col = c(colors[c(12,9,10,11)], "gray40"), pch = c(19,19,19,19,1),
       cex= 0.8, bty = "n")


dev.off()
## for text:

### Significant with robust version #### 
table(m_robust$omics_id[!unknown], m_robust$q_value[!unknown] < 0.05)
# 
#   FALSE  TRUE
# 5 10726  2537
# 1   358   159
# 2   483   163
# 3    55    11
# 4    36     3

table(m_robust$q_value[unknown] < 0.05)

# FALSE  TRUE 
# 2253   503 


##### Gene set enrichment ######

## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# pull
metadata_GO_bp <- dbGetQuery(con, "SELECT metadata_type, metadata_value, metadata.biomolecule_id, standardized_name 
                              FROM metadata 
                              INNER JOIN biomolecules ON metadata.biomolecule_id = biomolecules.biomolecule_id
                             WHERE metadata_type = 'GO_biological_process'")
metadata_class <- dbGetQuery(con, "SELECT metadata_type, metadata_value, metadata.biomolecule_id, standardized_name 
                              FROM metadata 
                              INNER JOIN biomolecules ON metadata.biomolecule_id = biomolecules.biomolecule_id
                             WHERE metadata_type = 'Larger_class' OR metadata_type = 'Lipid Class'")

metadata_geneName <- dbGetQuery(con, "SELECT metadata_type, metadata_value, metadata.biomolecule_id, standardized_name 
                              FROM metadata 
                              INNER JOIN biomolecules ON metadata.biomolecule_id = biomolecules.biomolecule_id
                             WHERE metadata_type = 'gene_name'")

# disconnect
dbDisconnect(con) 

# 

metadata_class <- metadata_class[!(metadata_class$metadata_value == ""),]
metadata_class$metadata_value <- paste(metadata_class$metadata_value,"[Metabolite Class]")


##### make reference set with go term bp and class info #### 
metadata_GO_bp_class<- rbind(metadata_GO_bp, metadata_class)


class_and_GO_bp <- make_reference_sets(metadata_GO_bp_class$metadata_value, metadata_GO_bp_class$biomolecule_id)
class_and_GO_bp_standardized_names <- make_reference_sets(metadata_GO_bp_class$metadata_value, metadata_GO_bp_class$standardized_name)

#### performing enrichment based Robust LR test significant for COVID #### 
enrichment_significant_up_COVID_robust <- enrichment(m_robust$biomolecule_id[m_robust$q_value < 0.05 & m_robust$FC > 1], class_and_GO_bp, m_robust$biomolecule_id)
enrichment_significant_down_COVID_robust <- enrichment(m_robust$biomolecule_id[m_robust$q_value < 0.05 & m_robust$FC < -1], class_and_GO_bp, m_robust$biomolecule_id)

length(m_robust$biomolecule_id[m_robust$q_value < 0.05 & m_robust$FC > 1]) #594
length(m_robust$biomolecule_id[m_robust$q_value < 0.05 & m_robust$FC < -1]) #384

enrichment_significant_up_COVID_robust_std <- enrichment(m_robust$standardized_name[m_robust$q_value < 0.05 & m_robust$FC > 1], class_and_GO_bp_standardized_names, m_robust$standardized_name)
enrichment_significant_down_COVID_robust_std <- enrichment(m_robust$standardized_name[m_robust$q_value < 0.05 & m_robust$FC < -1], class_and_GO_bp_standardized_names, m_robust$standardized_name)

write.csv(cbind(enrichment_significant_up_COVID_robust_std, enrichment_significant_down_COVID_robust_std), "data/enrichment_significant_COVID_robustlm.csv")

enrichment_significant_up_COVID_protein <- enrichment(m_robust$standardized_name[m_robust$q_value < 0.05 & m_robust$FC > 0 & m_robust$omics_id == 1], class_and_GO_bp_standardized_names, m_robust$standardized_name)
enrichment_significant_down_COVID_protein <- enrichment(m_robust$standardized_name[m_robust$q_value < 0.05 & m_robust$FC < 0 & m_robust$omics_id == 1], class_and_GO_bp_standardized_names, m_robust$standardized_name)

enrichment_significant_up_COVID_transcipts <- enrichment(m_robust$standardized_name[m_robust$q_value < 0.05 & m_robust$FC > 0 & m_robust$omics_id == 5], class_and_GO_bp_standardized_names, m_robust$standardized_name)
enrichment_significant_down_COVID_transcripts <- enrichment(m_robust$standardized_name[m_robust$q_value < 0.05 & m_robust$FC < 0 & m_robust$omics_id == 5], class_and_GO_bp_standardized_names, m_robust$standardized_name)

write.csv(cbind(enrichment_significant_up_COVID_protein, enrichment_significant_down_COVID_protein, enrichment_significant_up_COVID_transcipts, enrichment_significant_down_COVID_transcripts), file = "data/differences_COVID_enrichment_proteins_transcipts_standardized_names_robustlm.csv")

pdf("plots/Fig2_KAO_GO_BP_Enrichment_based_on_UP_COVID_v2_robustlm.pdf", width = 7, height = 6)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)

enrichment_significant_down_COVID_robust [order(enrichment_significant_down_COVID_robust [,2]),][1:7,]

barplot(-log10(enrichment_significant_up_COVID_robust [order(enrichment_significant_up_COVID_robust [,2]),2][1:7][7:1]), 
        xlab ="-log10(q-value)", 
        cex.names = 0.9, horiz = T, col = 1, 
        names = enrichment_significant_up_COVID_robust [order(enrichment_significant_up_COVID_robust [,2]),1][1:7][7:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significantly up with COVID")
dev.off()

pdf("plots/Fig2_KAO_GO_BP_Enrichment_based_on_DOWN_COVID_v2_robustlm.pdf", width = 7, height = 2.5)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)

barplot(-log10(enrichment_significant_down_COVID_robust[order(enrichment_significant_down_COVID_robust[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 3, 
        names = enrichment_significant_down_COVID_robust[order(enrichment_significant_down_COVID_robust[,2]),1][1:5][5:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significantly down with COVID")

dev.off()



pdf("plots/enriched_in_covid_volcano_plots_3_robust.pdf", useDingbats = F, height = 7, width = 5)
par(mfrow = c(3,2))
for(i in c(11, 12, 22,23,133)){
  up <- class_and_GO_bp[[enrichment_significant_up_COVID_robust[order(enrichment_significant_up_COVID_robust$pvalue),][i,1]]]
  
  plot(m_robust$FC[!unknown], -log(m_robust$q_value)[!unknown], 
       pch = 19, col = "gray70",
       ylim = c(0, 40), xlim = c(-4, 6),
       main = paste(enrichment_significant_up_COVID_robust[order(enrichment_significant_up_COVID_robust$pvalue),][i,1]), 
       ylab = "-log(adjusted p-value of likelyhood ratio)",
       xlab = "log2(abundance COVID-19 positive/COVID-19 negative)",
       las = 1, bty = "l", cex.axis = 1.2)
  points(m_robust$FC[m_robust$biomolecule_id %in% up], -log(m_robust$q_value)[m_robust$biomolecule_id %in% up], 
         pch = 19, col = colors[c(12,9,10,11,11)][m_robust$omics_id[m_robust$biomolecule_id %in% up]])
  
  
}


for(i in c(1)){
  down <- class_and_GO_bp[[enrichment_significant_down_COVID_robust[order(enrichment_significant_down_COVID_robust$pvalue),][i,1]]]
  
  plot(m_robust$FC[!unknown], -log(m_robust$q_value)[!unknown], 
       pch = 19, col = "gray70",
       ylim = c(0, 40), xlim = c(-4, 6),
       main = paste(enrichment_significant_down_COVID_robust[order(enrichment_significant_down_COVID_robust$pvalue),][i,1]), 
       ylab = "-log(adjusted p-value of likelyhood ratio)",
       xlab = "log2(abundance COVID-19 positive/COVID-19 negative)",
       las = 1, bty = "l", cex.axis = 1.2)
  points(m_robust$FC[m_robust$biomolecule_id %in% down], -log(m_robust$q_value)[m_robust$biomolecule_id %in% down], 
         pch = 19, col = colors[c(12,9,10,11,11)][m_robust$omics_id[m_robust$biomolecule_id %in% down]])
  
}

dev.off()

#### HFDs significance ##### 

# connect 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# pull pvalues 

pvalues_HFD_robust <- dbGetQuery(con, "SELECT * FROM pvalues 
          INNER JOIN biomolecules on biomolecules.biomolecule_id = pvalues.biomolecule_id
          WHERE pvalues.formula = 'Hospital_free_days_45 ~ normalized_abundance + Age_less_than_90 + Gender vs. Hospital_free_days_45 ~ Age_less_than_90 + Gender' AND pvalues.test = 'LR_test_robust' AND biomolecules.keep = 1
                  ")
# disconnect
dbDisconnect(con) 

## combine with other pvalues table robust 

pvalues_HFD_robust <- pvalues_HFD_robust[,-grep('biomolecule_id', names(pvalues_HFD_robust))[2]]

m2_robust <- merge(m_robust, pvalues_HFD_robust, by = 'biomolecule_id')

unknown_m2_robust <- grepl("nknown", m2_robust$standardized_name.x)

table(m2_robust$omics_id.x[!unknown_m2_robust], m2_robust$q_value.y[!unknown_m2_robust] < 0.05)

#    FALSE TRUE
# 5  6980 6283
# 1   322  195
# 2   423  223
# 3    42   24
# 4    27   12

## direction of change for HFD ## 
UP_with_severity <- colMeans(df_wide_exprs[bin_values <3, ]) >  colMeans(df_wide_exprs[bin_values > 7, ])

direction_HFD <- data.frame(biomolecule_id = sub("normalized_abundance.", "", names(df_wide_exprs)), UP_with_severity = UP_with_severity)

m3_robust <- merge(m2_robust, direction_HFD, by = "biomolecule_id")
m4_robust <- merge(m3_robust, metadata_geneName, by = "biomolecule_id", all.x = T)


#### enrichment analysis for p-values of HFD ####
enrichment_significant_up_HFD_robust <- enrichment(m3_robust$biomolecule_id[m3_robust$q_value.y < 0.05 & m3_robust$UP_with_severity], class_and_GO_bp, m3_robust$biomolecule_id)
enrichment_significant_down_HFD_robust <- enrichment(m3_robust$biomolecule_id[m3_robust$q_value.y < 0.05 & !m3_robust$UP_with_severity], class_and_GO_bp, m3_robust$biomolecule_id)

pdf("plots/Fig2_KAO_GO_BP_Enrichment_based_on_UP_Severity_v2_robust.pdf", width = 7, height = 6)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)

barplot(-log10(enrichment_significant_up_HFD_robust [order(enrichment_significant_up_HFD_robust [,2]),2][1:30][30:1]), 
        xlab ="-log10(q-value)", 
        cex.names = 0.9, horiz = T, col = "black", 
        names = enrichment_significant_up_HFD_robust [order(enrichment_significant_up_HFD_robust [,2]),1][1:30][30:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significantly up with severity")
dev.off()

pdf("plots/Fig2_KAO_GO_BP_Enrichment_based_on_DOWN_severity_v2_robust.pdf", width = 7, height = 6)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)

barplot(-log10(enrichment_significant_down_HFD_robust[order(enrichment_significant_down_HFD_robust[,2]),2][1:30][30:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = "gray80", 
        names = enrichment_significant_down_HFD_robust[order(enrichment_significant_down_HFD_robust[,2]),1][1:30][30:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significantly down with severity")

dev.off()


#### BOTH HFD and COVID #### 
table(m2_robust$q_value.x[!unknown_m2_robust] < 0.05, m2_robust$q_value.y[!unknown_m2_robust] < 0.05, m2_robust$omics_id.x[!unknown_m2_robust])
# , ,  = 5
# 
# 
# FALSE TRUE
# FALSE  5291 5435
# TRUE   1689  848
# 
# , ,  = 1
# 
# 
# FALSE TRUE
# FALSE   264   94
# TRUE     58  101
# 
# , ,  = 2
# 
# 
# FALSE TRUE
# FALSE   341  142
# TRUE     82   81
# 
# , ,  = 3
# 
# 
# FALSE TRUE
# FALSE    35   20
# TRUE      7    4
# 
# , ,  = 4
# 
# 
# FALSE TRUE
# FALSE    25   11
# TRUE      2    1


#### Venn diagram both #### 

library(VennDiagram)

#### Heatmap with HFD and COVID #####

both_robust <- m2_robust$biomolecule_id[!grepl("nknown", m3_robust$standardized_name.x) & (m3_robust$q_value.x < 0.05 & m3_robust$FC > 1 & m3_robust$q_value.y < 0.05 & m3_robust$UP_with_severity) | (m3_robust$q_value.x < 0.05 & m3_robust$FC < -1 & m3_robust$q_value.y < 0.05 & !m3_robust$UP_with_severity)] 

names(df_wide_exprs) <- sub("normalized_abundance.", "", names(df_wide_exprs)) 
t_df_wide_exprs <- t(df_wide_exprs)

table( row.names(t_df_wide_exprs) %in% both_robust)
# FALSE  TRUE 
# 16806   454

table(m3_robust$omics_id.x, !grepl("nknown", m3_robust$standardized_name.x) & (m3_robust$q_value.x < 0.05 & m3_robust$FC > 1 & m3_robust$q_value.y < 0.05 & m3_robust$UP_with_severity) | (m3_robust$q_value.x < 0.05 & m3_robust$FC < -1 & m3_robust$q_value.y < 0.05 & !m3_robust$UP_with_severity))
# FALSE  TRUE
# 5 12876   387
# 1   495    22
# 2  3312    45
# 3   111     0
# 4    39     0


#### Elastic Net results #### 

EN <- do.call(rbind, lapply(list.files("./data/ElasticNet", "tsv", recursive = T, full.names=TRUE), read.delim, stringsAsFactors = F))
names(EN) <- c("standardized_name", "EN_coefficient")

EN_merge_robust <- merge(m4_robust, EN, by.x = "standardized_name.x", by.y = 1, all.x = T)

EN_merge_robust <- EN_merge_robust[,c(1,2,8,9,10,18,24,25,29,31,33)]
names(EN_merge_robust) <- sub(".x", ".COVID", names(EN_merge_robust), fixed = T)
names(EN_merge_robust) <- sub(".y", ".HFD", names(EN_merge_robust), fixed = T)

table(!is.na(EN_merge_robust$EN_coefficient) & !(grepl("nknown", EN_merge_robust$standardized_name.COVID)), EN_merge_robust$omics_id.COVID)

## Joint Venn Diagram 

##### Common with robust methods #### 
biomolecules_COVID_robust <- EN_merge_robust$biomolecule_id[EN_merge_robust$q_value.COVID < 0.05]
biomolecules_HFD_robust <- EN_merge_robust$biomolecule_id[EN_merge_robust$q_value.HFD < 0.05]
biomolecule_EN_robust <- EN_merge_robust$biomolecule_id[!is.na(EN_merge_robust$EN_coefficient)]

biomolecule_intersect_robust <- Reduce(intersect, list(biomolecule_EN_robust, biomolecules_COVID_robust, biomolecules_HFD_robust))

length(biomolecule_intersect_robust) #255

## with same directions 

#up
biomolecules_COVID_up_robust <- EN_merge_robust$biomolecule_id[EN_merge_robust$q_value.COVID < 0.05 & EN_merge_robust$FC > 0]
biomolecules_HFD_up_robust <- EN_merge_robust$biomolecule_id[EN_merge_robust$q_value.HFD < 0.05 & EN_merge_robust$UP_with_severity]
biomolecule_EN_up_robust <- EN_merge_robust$biomolecule_id[!is.na(EN_merge_robust$EN_coefficient) & EN_merge_robust$EN_coefficient < 0]

biomolecule_intersect_up_robust <- Reduce(intersect, list(biomolecule_EN_up_robust, biomolecules_COVID_up_robust, biomolecules_HFD_up_robust))

length(biomolecule_intersect_up_robust) #142

#down
biomolecules_COVID_down_robust <- EN_merge_robust$biomolecule_id[EN_merge_robust$q_value.COVID < 0.05 & EN_merge_robust$FC < 0]
biomolecules_HFD_down_robust <- EN_merge_robust$biomolecule_id[EN_merge_robust$q_value.HFD < 0.05 & !EN_merge_robust$UP_with_severity]
biomolecule_EN_down_robust <- EN_merge_robust$biomolecule_id[!is.na(EN_merge_robust$EN_coefficient) & EN_merge_robust$EN_coefficient > 0]

biomolecule_intersect_down_robust <- Reduce(intersect, list(biomolecule_EN_down_robust, biomolecules_COVID_down_robust, biomolecules_HFD_down_robust))

length(biomolecule_intersect_down_robust) #85

biomolecule_interest_directional_robust <- c(biomolecule_intersect_down_robust, biomolecule_intersect_up_robust)
biomolecule_unkowns_robust <- EN_merge_robust$biomolecule_id[grepl("nknown", EN_merge_robust$standardized_name.COVID)]

length(setdiff(biomolecule_interest_directional_robust, biomolecule_unkowns_robust)) #140

pdf("plots/overlap_COVID_HFD_EN_robust.pdf")

venn.plot <- draw.triple.venn(
  area1 = length(c(biomolecules_COVID_down_robust, biomolecules_COVID_up_robust)),
  area2 = length(c(biomolecules_HFD_down_robust, biomolecules_HFD_up_robust)),
  area3 = length(c(biomolecule_EN_down_robust, biomolecule_EN_up_robust)), 
  n12 = length(c(intersect(biomolecules_COVID_down_robust,biomolecules_HFD_down_robust), intersect(biomolecules_COVID_up_robust, biomolecules_HFD_up_robust))), 
  n23 = length(c(intersect(biomolecules_HFD_down_robust, biomolecule_EN_down_robust), intersect(biomolecules_HFD_up_robust, biomolecule_EN_up_robust))),
  n13 = 
    length(c(intersect(biomolecules_COVID_down_robust, biomolecule_EN_down_robust), intersect(biomolecules_COVID_up_robust, biomolecule_EN_up_robust))),
  n123 = length(biomolecule_interest_directional_robust),
  fontfamily = rep("",7),
  cat.fontfamily = rep("", 3),
  main.fontfamily = "",
  fill = c("gray40", "gray40", "gray40"),
  category= c( "Associated with COVID-19", "Associated with HFD-45", "Elastic Net feature selection"), 
  main = "High interest features")
dev.off()

#### EN heatmap #####

annotation_row <- data.frame(omics_id = m2_robust$omics_id.x)
row.names(annotation_row) <- m2_robust$biomolecule_id

annotation_col <- data.frame(COVID = df_wide_all$COVID, 
                             HFD = bin_values)
row.names(annotation_col) <- row.names(df_wide_all)

# annotation colors 
annotation_colors <- list(omics_id = colors[c(12,9,10,11,11)], 
                          COVID = colors[c(3,1)], 
                          HFD = colorRampPalette(c("black", "gray75"))(9))

names(annotation_colors[["omics_id"]]) <- as.character(levels(annotation_row$omics_id))
names(annotation_colors[["COVID"]]) <- as.character(levels(annotation_col$COVID))
names(annotation_colors[["HFD"]]) <- as.character(levels(annotation_col$HFD))

pdf("plots/heatmap_bothHFD_COVID_significant_Fig2_no_unknowns_includes_ENFilter_robust.pdf")
pheatmap((t_df_wide_exprs[row.names(t_df_wide_exprs) %in% setdiff(biomolecule_interest_directional_robust, biomolecule_unkowns_robust), ]), 
         scale = "row",
         breaks = c(-10, -5, seq(-2,2,length.out = 16),5, 10),
         color = colorRampPalette(c(3,"gray90",1))(20),
         clustering_method = "ward.D2",
         annotation_col = annotation_col, annotation_row =annotation_row,
         annotation_colors = annotation_colors,
         show_rownames = F,
         show_colnames = F
)
dev.off()

### for robust heatmap
annotation_row <- data.frame(omics_id = m2_robust$omics_id.x)
row.names(annotation_row) <- m2_robust$biomolecule_id

annotation_col <- data.frame(COVID = df_wide_all$COVID, 
                             HFD = bin_values)
row.names(annotation_col) <- row.names(df_wide_all)

# annotation colors 
annotation_colors <- list(omics_id = colors[c(12,9,10,11,11)], 
                          COVID = colors[c(3,1)], 
                          HFD = colorRampPalette(c("black", "gray75"))(9))

names(annotation_colors[["omics_id"]]) <- as.character(levels(annotation_row$omics_id))
names(annotation_colors[["COVID"]]) <- as.character(levels(annotation_col$COVID))
names(annotation_colors[["HFD"]]) <- as.character(levels(annotation_col$HFD))

pdf("plots/heatmap_bothHFD_COVID_significant_Fig2_no_unknowns_includes_ENFilter_robust.pdf")
pheatmap((t_df_wide_exprs[row.names(t_df_wide_exprs) %in% setdiff(biomolecule_interest_directional_robust, biomolecule_unkowns_robust), ]), 
         scale = "row",
         breaks = c(-10, -5, seq(-2,2,length.out = 16),5, 10),
         color = colorRampPalette(c(3,"gray90",1))(20),
         clustering_method = "ward.D2",
         annotation_col = annotation_col, annotation_row =annotation_row,
         annotation_colors = annotation_colors,
         show_rownames = F,
         show_colnames = F
)
dev.off()


### enrichment 
enrichment_significant_both_robust <- enrichment(biomolecule_interest_directional_robust, class_and_GO_bp, m3_robust$biomolecule_id)
enrichment_intersect_down_robust <- enrichment(biomolecule_intersect_down_robust, class_and_GO_bp, m3_robust$biomolecule_id)
enrichment_intersect_up_robust <- enrichment(biomolecule_intersect_up_robust, class_and_GO_bp, m3_robust$biomolecule_id)

enrichment_intersect_up_robust[5,]

enrichment_intersect_down_std_robust <- enrichment(EN_merge_robust$standardized_name.COVID[EN_merge_robust$biomolecule_id %in% biomolecule_intersect_down_robust], class_and_GO_bp_standardized_names, EN_merge_robust$standardized_name.COVID)
enrichment_intersect_up_std_robust <- enrichment(EN_merge_robust$standardized_name.COVID[EN_merge_robust$biomolecule_id %in% biomolecule_intersect_up_robust], class_and_GO_bp_standardized_names, EN_merge_robust$standardized_name.COVID)


pdf("plots/Fig2_KAO_GO_BP_Enrichment_based_on_UP_intersect_robust.pdf", width = 7, height =5.5)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
par(mfrow = c(2,1))
barplot(-log10(enrichment_intersect_up_robust [order(enrichment_intersect_up_robust[,2]),2][1:7][7:1]), 
        xlab ="-log10(q-value)", 
        cex.names = 0.9, horiz = T, col = 1, 
        names = enrichment_intersect_up_robust [order(enrichment_intersect_up_robust [,2]),1][1:7][7:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significantly up with COVID status and severity")

barplot(-log10(enrichment_intersect_down_robust[order(enrichment_intersect_down_robust[,2]),2][1:7][7:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 3, 
        names = enrichment_intersect_down_robust[order(enrichment_intersect_down_robust[,2]),1][1:7][7:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significantly down with COVID status and severity")

dev.off()

write.csv(EN_merge_robust[EN_merge_robust$biomolecule_id %in% biomolecule_interest_directional_robust,], "intersect_features_robust.csv")
# 
# df_forScott<- t_df_wide_exprs[row.names(t_df_wide_exprs) %in% setdiff(biomolecule_interest_directional, biomolecule_unkowns), ]
# colnames(df_forScott) <- df_wide_all$sample_id
# write.csv(df_forScott, "data/Fig2d_feature_table.csv")

write.csv(cbind(enrichment_intersect_up_std_robust, enrichment_intersect_down_std_robust), "data/enrichment_intersect_COVID_and_sevity_robust.csv")
