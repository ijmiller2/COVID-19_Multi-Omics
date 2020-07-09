##### 14_KAO_figure_2_version_1.R ##### 

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


df_metabolites<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, metabolomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45
           FROM metabolomics_measurements
           INNER JOIN metabolomics_runs ON metabolomics_runs.replicate_id = metabolomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = metabolomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = metabolomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")

df_lipids<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, lipidomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45
           FROM lipidomics_measurements
           INNER JOIN lipidomics_runs ON lipidomics_runs.replicate_id = lipidomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = lipidomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = lipidomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_proteins<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, proteomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45
           FROM proteomics_measurements
           INNER JOIN proteomics_runs ON proteomics_runs.replicate_id = proteomics_measurements.replicate_id
           INNER JOIN rawfiles ON rawfiles.rawfile_id = proteomics_runs.rawfile_id
           INNER JOIN deidentified_patient_metadata ON deidentified_patient_metadata.sample_id = rawfiles.sample_id
           INNER JOIN biomolecules on biomolecules.biomolecule_id = proteomics_measurements.biomolecule_id
           WHERE rawfiles.keep = 1  
           AND biomolecules.keep = '1'
           ")


df_transcripts<- dbGetQuery(con, "SELECT deidentified_patient_metadata.sample_id, normalized_abundance, transcriptomics_measurements.biomolecule_id, COVID, ICU_1, Hospital_free_days_45
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

df_wide_exprs <- df_wide_all[,-c(1:4)]


##### Priciple component analysis #######
pca1 <- prcomp(df_wide_exprs)

## Creating color scale for HFD
breaks = seq(0,45, by =5)
bin_values <- .bincode(df_wide_all$Hospital_free_days_45, breaks, include.lowest = T)
col_HFD<- colorRampPalette(c("black", "gray75"))(9)[bin_values]

## creating ellipse with 95% CI for top and bottom bin values
lvs2 <- levels(as.factor(bin_values))
pts.array2 <- array(0, dim = c(100, 2, length(lvs2)))

for(i in c(1,9)) {
  inx <- bin_values == lvs2[i]
  groupVar <- var(cbind(pca1$x[inx, 1], pca1$x[inx,2]), na.rm = T )
  groupMean <- cbind(mean(pca1$x[inx, 1], na.rm= T), mean(pca1$x[inx, 2], na.rm =T))
  pts.array2[,,i] <- ellipse(groupVar, centre = groupMean, level = 0.95, npoints = 100)
}


## Plotting PCA 
pdf("plots/Fig2_PCA_all_omes_HFD_color.pdf", width = 5, height = 5, useDingbats = F)
plt_original <- par()$plt
par(mgp = c(1.8, 0.7, 0), tcl = -0.3)
par("plt" = c(0.92,0.97, 0.22, 0.43))
plot(breaks, type ="n", bty = "n", xaxt = "n", xlab = "", ylab = "HFD-45", las =1 , )
for(i in 1:length(breaks[-1])){
  rect(0,breaks[i],10,breaks[i+1], col = colorRampPalette(c("black", "gray75"))(9)[i], border =NA)
}

par(new = T)
par("plt" = plt_original)

plot(pca1$x, col= col_HFD, pch =c(17,19)[df_wide_all$COVID+1], cex = 1.5, lwd = 3, bty = "l",
     ylim = c(-125, 125), xlim = c(-175, 150), las =1 ,
     main = "Priciple component anlaysis - all omes",
     xlab = paste("PC1 (", round(summary(pca1)$importance[2,1]*100, digits = 2), "%)", sep = ""),
     ylab = paste("PC2 (", round(summary(pca1)$importance[2,2]*100, digits = 2), "%)", sep = ""))


polygon(pts.array2[,,1], border = "black", lty = 2, lwd = 2)
polygon(pts.array2[,,9], border = "gray75", lty = 2, lwd = 2)

legend("bottomleft", pch = c(19, 17), c("COVID-19", "non-COVID-19"), bty = "n")


dev.off()

     
### Fig 2B - volcano plot COVID-status #### 

# connect 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# pull pvalues 
pvalues_COVID <- dbGetQuery(con, "SELECT * FROM pvalues 
          INNER JOIN biomolecules on biomolecules.biomolecule_id = pvalues.biomolecule_id
          WHERE pvalues.comparison = 'COVID_vs_NONCOVID' AND biomolecules.keep = 1
                  ")

# disconnect
dbDisconnect(con) 

## load in log2 fold change values 

FC <- read.csv("data/COVID_fc_by_biomolecule_ID.csv")

## load in transcript pvalues 

transcript_pvalues <- read.delim("P:/All_20200428_COVID_plasma_multiomics/Transcriptomics/AllGenes.txt", sep = " ")
transcript_pvalues <- data.frame(standardized_name = row.names(transcript_pvalues), PP = 1.000000001-transcript_pvalues[,1], stringsAsFactors = F)

## merging pvalues from database with DESeq results

pvalues_COVID <- pvalues_COVID[,-grep('biomolecule_id', names(pvalues_COVID))[2]]
pvalues_COVID_deseq <- merge(pvalues_COVID, transcript_pvalues, by = "standardized_name", all.x = T)

## replacing original q-values based on log-likelyhood with DEseq results (1-PP)

pvalues_COVID_deseq$q_value[!is.na(pvalues_COVID_deseq$PP)] <- pvalues_COVID_deseq$PP[!is.na(pvalues_COVID_deseq$PP)] 

## merge pvalues and FC by biomolecule_id

m <- merge(pvalues_COVID_deseq, FC, by.x = 'biomolecule_id', by.y = 'biomolecule_id')
m$omics_id <- as.factor(m$omics_id)
m$omics_id <- relevel(m$omics_id, ref = 5)
m <- m[order(m$omics_id), ]

unknown <- grepl("nknown", m$standardized_name)
table(unknown)

## Plot Volcano plot 
pdf("plots/KAO_fig_2A_Volcano_plot_v2.pdf", useDingbats = F, width = 5, height = 5)

par(mgp = c(1.8, 0.7, 0), tcl = -0.3)
plot(m$FC[unknown], -log(m$q_value)[unknown],
     pch = 1, col = "gray40",
     ylim = c(0, 34), xlim = c(-4, 6),
     main = "Effect of COVID", 
     ylab = "-log(adjusted p-value)",
     xlab = "log2(COVID-19/non-COVID-19)",
     las = 1, bty = "l", cex.axis = 1.2)
points(m$FC[!unknown], -log(m$q_value)[!unknown], 
     pch = 19, col = colors[c(12,9,10,11,11)][m$omics_id[!unknown]])
     
legend("topleft", c("transcripts", "proteins", "lipids", "small molecules", "unidentified features"), 
       col = c(colors[c(12,9,10,11)], "gray40"), pch = c(19,19,19,19,1),
       cex= 0.8, bty = "n")

dev.off()

## for text:

table(m$omics_id[!unknown],m$q_value[!unknown] < 0.05)

# FALSE  TRUE
# 5 10726  2537
# 1   371   146
# 2   478   168
# 3    58     8
# 4    34     5

table(m$q_value[unknown] < 0.05)

# FALSE  TRUE 
# 2245   511 

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

#### performing enrichment based LR test significant for COVID #### 
enrichment_significant_up_COVID <- enrichment(m$biomolecule_id[m$q_value < 0.05 & m$FC > 1], class_and_GO_bp, m$biomolecule_id)
enrichment_significant_down_COVID <- enrichment(m$biomolecule_id[m$q_value < 0.05 & m$FC < -1], class_and_GO_bp, m$biomolecule_id)

length(m$biomolecule_id[m$q_value < 0.05 & m$FC > 1]) #586
length(m$biomolecule_id[m$q_value < 0.05 & m$FC < -1]) #383

enrichment_significant_up_COVID_std <- enrichment(m$standardized_name[m$q_value < 0.05 & m$FC > 1], class_and_GO_bp_standardized_names, m$standardized_name)
enrichment_significant_down_COVID_std <- enrichment(m$standardized_name[m$q_value < 0.05 & m$FC < -1], class_and_GO_bp_standardized_names, m$standardized_name)


write.csv(cbind(enrichment_significant_up_COVID_std, enrichment_significant_down_COVID_std), "data/enrichment_significant_COVID.csv")

enrichment_significant_up_COVID_protein <- enrichment(m$standardized_name[m$q_value < 0.05 & m$FC > 0 & m$omics_id == 1], class_and_GO_bp_standardized_names, m$standardized_name)
enrichment_significant_down_COVID_protein <- enrichment(m$standardized_name[m$q_value < 0.05 & m$FC < 0 & m$omics_id == 1], class_and_GO_bp_standardized_names, m$standardized_name)

enrichment_significant_up_COVID_transcipts <- enrichment(m$standardized_name[m$q_value < 0.05 & m$FC > 0 & m$omics_id == 5], class_and_GO_bp_standardized_names, m$standardized_name)
enrichment_significant_down_COVID_transcripts <- enrichment(m$standardized_name[m$q_value < 0.05 & m$FC < 0 & m$omics_id == 5], class_and_GO_bp_standardized_names, m$standardized_name)

write.csv(cbind(enrichment_significant_up_COVID_protein, enrichment_significant_down_COVID_protein, enrichment_significant_up_COVID_transcipts, enrichment_significant_down_COVID_transcripts), file = "data/differences_COVID_enrichment_proteins_transcipts_standardized_names.csv")

pdf("plots/Fig2_KAO_GO_BP_Enrichment_based_on_UP_COVID_v2.pdf", width = 7, height = 6)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)

enrichment_significant_down_COVID [order(enrichment_significant_down_COVID [,2]),][1:30,]

barplot(-log10(enrichment_significant_up_COVID [order(enrichment_significant_up_COVID [,2]),2][1:30][30:1]), 
        xlab ="-log10(q-value)", 
        cex.names = 0.9, horiz = T, col = 1, 
        names = enrichment_significant_up_COVID [order(enrichment_significant_up_COVID [,2]),1][1:30][30:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significantly up with COVID")
dev.off()

pdf("plots/Fig2_KAO_GO_BP_Enrichment_based_on_DOWN_COVID_v2.pdf", width = 7, height = 2.5)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)

barplot(-log10(enrichment_significant_down_COVID[order(enrichment_significant_down_COVID[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 3, 
        names = enrichment_significant_down_COVID[order(enrichment_significant_down_COVID[,2]),1][1:5][5:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significantly down with COVID")

dev.off()



pdf("plots/enriched_in_covid_volcano_plots_3.pdf", useDingbats = F, height = 7, width = 5)
par(mfrow = c(3,2))
for(i in c(10, 11, 21,22,31)){
up <- class_and_GO_bp[[enrichment_significant_up_COVID[order(enrichment_significant_up_COVID$pvalue),][i,1]]]

plot(m$FC[!unknown], -log(m$q_value)[!unknown], 
     pch = 19, col = "gray70",
     ylim = c(0, 40), xlim = c(-4, 6),
     main = paste(enrichment_significant_up_COVID[order(enrichment_significant_up_COVID$pvalue),][i,1]), 
     ylab = "-log(adjusted p-value of likelyhood ratio)",
     xlab = "log2(abundance COVID-19 positive/COVID-19 negative)",
     las = 1, bty = "l", cex.axis = 1.2)
points(m$FC[m$biomolecule_id %in% up], -log(m$q_value)[m$biomolecule_id %in% up], 
       pch = 19, col = colors[c(12,9,10,11,11)][m$omics_id[m$biomolecule_id %in% up]])


}


for(i in c(1)){
  down <- class_and_GO_bp[[enrichment_significant_down_COVID[order(enrichment_significant_down_COVID$pvalue),][i,1]]]
  
  plot(m$FC[!unknown], -log(m$q_value)[!unknown], 
       pch = 19, col = "gray70",
       ylim = c(0, 40), xlim = c(-4, 6),
       main = paste(enrichment_significant_down_COVID[order(enrichment_significant_down_COVID$pvalue),][i,1]), 
       ylab = "-log(adjusted p-value of likelyhood ratio)",
       xlab = "log2(abundance COVID-19 positive/COVID-19 negative)",
       las = 1, bty = "l", cex.axis = 1.2)
  points(m$FC[m$biomolecule_id %in% down], -log(m$q_value)[m$biomolecule_id %in% down], 
         pch = 19, col = colors[c(12,9,10,11,11)][m$omics_id[m$biomolecule_id %in% down]])
  
}

dev.off()

#### HFDs significance ##### 

# connect 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# pull pvalues 
pvalues_HFD <- dbGetQuery(con, "SELECT * FROM pvalues 
          INNER JOIN biomolecules on biomolecules.biomolecule_id = pvalues.biomolecule_id
          WHERE pvalues.formula = 'Hospital_free_days_45 ~ normalized_abundance + Age_less_than_90 + Gender vs. Hospital_free_days_45 ~ Age_less_than_90 + Gender' AND biomolecules.keep = 1
                  ")

# disconnect
dbDisconnect(con) 

## combine with other pvalues table 
pvalues_HFD <- pvalues_HFD[,-grep('biomolecule_id', names(pvalues_HFD))[2]]

m2 <- merge(m, pvalues_HFD, by = 'biomolecule_id')

unknown_m2 <- grepl("nknown", m2$standardized_name.x)
## for text:
table(m2$omics_id.x[!unknown_m2], m2$q_value.y[!unknown_m2] < 0.05)

# FALSE TRUE
# 1   328  189
# 2   428  218
# 3    43   23
# 4    27   12
# 5  7061 6202

table(pvalues_HFD$q_value[unknown_HFD] < 0.05)

# FALSE  TRUE 
# 1992   764

## direction of change for HFD ## 
UP_with_severity <- colMeans(df_wide_exprs[bin_values <3, ]) >  colMeans(df_wide_exprs[bin_values > 7, ])

direction_HFD <- data.frame(biomolecule_id = sub("normalized_abundance.", "", names(df_wide_exprs)), UP_with_severity = UP_with_severity)

m3 <- merge(m2, direction_HFD, by = "biomolecule_id")
m4 <- merge(m3, metadata_geneName, by = "biomolecule_id", all.x = T)

#### enrichment analysis for p-values of HFD ####
enrichment_significant_up_HFD <- enrichment(m3$biomolecule_id[m3$q_value.y < 0.05 & m3$UP_with_severity], class_and_GO_bp, m3$biomolecule_id)
enrichment_significant_down_HFD <- enrichment(m3$biomolecule_id[m3$q_value.y < 0.05 & !m3$UP_with_severity], class_and_GO_bp, m3$biomolecule_id)

pdf("plots/Fig2_KAO_GO_BP_Enrichment_based_on_UP_Severity_v2.pdf", width = 7, height = 6)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)

barplot(-log10(enrichment_significant_up_HFD [order(enrichment_significant_up_HFD [,2]),2][1:30][30:1]), 
        xlab ="-log10(q-value)", 
        cex.names = 0.9, horiz = T, col = "black", 
        names = enrichment_significant_up_HFD [order(enrichment_significant_up_HFD [,2]),1][1:30][30:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significantly up with severity")
dev.off()

pdf("plots/Fig2_KAO_GO_BP_Enrichment_based_on_DOWN_severity_v2.pdf", width = 7, height = 6)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)

barplot(-log10(enrichment_significant_down_HFD[order(enrichment_significant_down_HFD[,2]),2][1:30][30:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = "gray80", 
        names = enrichment_significant_down_HFD[order(enrichment_significant_down_HFD[,2]),1][1:30][30:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significantly down with severity")

dev.off()


#### BOTH HFD and COVID #### 

table(m2$q_value.x[!unknown_m2] < 0.05, m2$q_value.y[!unknown_m2] < 0.05, m2$omics_id.x[!unknown_HFD])

# 
#  = 1
# FALSE TRUE
# FALSE   272   99
# TRUE     56   90
# 
#  = 2 
# FALSE TRUE
# FALSE   337  141
# TRUE     91   77
# 
#  = 3
# FALSE TRUE
# FALSE    37   21
# TRUE      6    2
# 
#  = 4
# FALSE TRUE
# FALSE    23   11
# TRUE      4    1
# 
#  = 5
# FALSE TRUE
# FALSE  5355 5371
# TRUE   1706  831

#### Venn diagram both #### 

library(VennDiagram)

m2$omics_id.x[m2$omics_id.x == 4] <-3


for (i in c(1,2,3,5)){

pdf(paste("plots/overlap_COVID_HFD_",i,".pdf", sep=""))
venn.plot <- draw.pairwise.venn(table(m2$q_value.y < 0.05 & m2$omics_id.x == i)[2], table(m2$q_value.x < 0.05 & m2$omics_id.x == i)[2], table(m2$q_value.x < 0.05 & m2$q_value.y < 0.05 & m2$omics_id.x == i)[2], fill = c("gray50",1),
                                category= c("Associated with HFD-45", "Associated with COVID-19"), main = "Overlap genus")
#grid.draw(venn.plot)
dev.off()

}

pdf("plots/overlap_COVID_HFD.pdf", width =4, height =4)
venn.plot <- draw.pairwise.venn(table(m2$q_value.y < 0.05)[2], table(m2$q_value.x < 0.05)[2], table(m2$q_value.x < 0.05 & m2$q_value.y < 0.05)[2], fill = c("gray50",1),
                                category= c("Associated with HFD-45", "Associated with COVID-19"), main = "Overlap genus")
dev.off()


#### Heatmap with HFD and COVID #####

both <- m2$biomolecule_id[!grepl("nknown", m3$standardized_name.x) & (m3$q_value.x < 0.05 & m3$FC > 1 & m3$q_value.y < 0.05 & m3$UP_with_severity) | (m3$q_value.x < 0.05 & m3$FC < -1 & m3$q_value.y < 0.05 & !m3$UP_with_severity)] 

names(df_wide_exprs) <- sub("normalized_abundance.", "", names(df_wide_exprs)) 
t_df_wide_exprs <- t(df_wide_exprs)

table( row.names(t_df_wide_exprs) %in% both)
# FALSE  TRUE 
# 16806   481 

table(m3$omics_id.x, !grepl("nknown", m3$standardized_name.x) & (m3$q_value.x < 0.05 & m3$FC > 1 & m3$q_value.y < 0.05 & m3$UP_with_severity) | (m3$q_value.x < 0.05 & m3$FC < -1 & m3$q_value.y < 0.05 & !m3$UP_with_severity))
# FALSE  TRUE
# 5 12885   378
# 1   496    21
# 2  3275    82
# 3   111     0
# 4    39     0


#### Elastic Net results #### 

EN <- do.call(rbind, lapply(list.files("./data/ElasticNet", "tsv", recursive = T, full.names=TRUE), read.delim, stringsAsFactors = F))
names(EN) <- c("standardized_name", "EN_coefficient")
EN_merge <- merge(m4, EN, by.x = "standardized_name.x", by.y = 1, all.x = T)

names(EN_merge)

EN_merge <- EN_merge[,c(1,2,8,9,10,18,24,25,29,31,33)]
names(EN_merge) <- sub(".x", ".COVID", names(EN_merge), fixed = T)
names(EN_merge) <- sub(".y", ".HFD", names(EN_merge), fixed = T)

table(!is.na(EN_merge$EN_coefficient) & !(grepl("nknown", EN_merge$standardized_name.COVID)), EN_merge$omics_id.COVID)

## Joint Venn Diagram 

biomolecules_COVID <- EN_merge$biomolecule_id[EN_merge$q_value.COVID < 0.05]
biomolecules_HFD <- EN_merge$biomolecule_id[EN_merge$q_value.HFD < 0.05]
biomolecule_EN <- EN_merge$biomolecule_id[!is.na(EN_merge$EN_coefficient)]

biomolecule_intersect <- Reduce(intersect, list(biomolecule_EN, biomolecules_COVID, biomolecules_HFD))

length(biomolecule_intersect) #248

## with same directions 

#up
biomolecules_COVID_up <- EN_merge$biomolecule_id[EN_merge$q_value.COVID < 0.05 & EN_merge$FC > 0]
biomolecules_HFD_up <- EN_merge$biomolecule_id[EN_merge$q_value.HFD < 0.05 & EN_merge$UP_with_severity]
biomolecule_EN_up <- EN_merge$biomolecule_id[!is.na(EN_merge$EN_coefficient) & EN_merge$EN_coefficient < 0]

biomolecule_intersect_up <- Reduce(intersect, list(biomolecule_EN_up, biomolecules_COVID_up, biomolecules_HFD_up))

length(biomolecule_intersect_up) #145

#down
biomolecules_COVID_down <- EN_merge$biomolecule_id[EN_merge$q_value.COVID < 0.05 & EN_merge$FC < 0]
biomolecules_HFD_down <- EN_merge$biomolecule_id[EN_merge$q_value.HFD < 0.05 & !EN_merge$UP_with_severity]
biomolecule_EN_down <- EN_merge$biomolecule_id[!is.na(EN_merge$EN_coefficient) & EN_merge$EN_coefficient > 0]

biomolecule_intersect_down <- Reduce(intersect, list(biomolecule_EN_down, biomolecules_COVID_down, biomolecules_HFD_down))

length(biomolecule_intersect_down) #74

biomolecule_interest_directional <- c(biomolecule_intersect_down, biomolecule_intersect_up)
biomolecule_unkowns <- EN_merge$biomolecule_id[grepl("nknown", EN_merge$standardized_name.COVID)]

length(setdiff(biomolecule_interest_directional, biomolecule_unkowns)) #132
pdf("plots/overlap_COVID_HFD_EN.pdf")

venn.plot <- draw.triple.venn(
  area1 = length(c(biomolecules_COVID_down, biomolecules_COVID_up)),
  area2 = length(c(biomolecules_HFD_down, biomolecules_HFD_up)),
  area3 = length(c(biomolecule_EN_down, biomolecule_EN_up)), 
  n12 = length(c(intersect(biomolecules_COVID_down,biomolecules_HFD_down), intersect(biomolecules_COVID_up, biomolecules_HFD_up))), 
  n23 = length(c(intersect(biomolecules_HFD_down, biomolecule_EN_down), intersect(biomolecules_HFD_up, biomolecule_EN_up))),
  n13 = 
    length(c(intersect(biomolecules_COVID_down, biomolecule_EN_down), intersect(biomolecules_COVID_up, biomolecule_EN_up))),
  n123 = length(biomolecule_interest_directional),
  fontfamily = rep("",7),
  cat.fontfamily = rep("", 3),
  main.fontfamily = "",
  fill = c("gray40", "gray40", "gray40"),
  category= c( "Associated with COVID-19", "Associated with HFD-45", "Elastic Net feature selection"), 
  main = "High interest features")
dev.off()

#### EN heatmap #####

annotation_row <- data.frame(omics_id = m2$omics_id.x)
row.names(annotation_row) <- m2$biomolecule_id

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

pdf("plots/heatmap_bothHFD_COVID_significant_Fig2_no_unknowns_includes_ENFilter.pdf")
pheatmap((t_df_wide_exprs[row.names(t_df_wide_exprs) %in% setdiff(biomolecule_interest_directional, biomolecule_unkowns), ]), 
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
enrichment_significant_both <- enrichment(biomolecule_interest_directional, class_and_GO_bp, m3$biomolecule_id)
enrichment_intersect_down <- enrichment(biomolecule_intersect_down, class_and_GO_bp, m3$biomolecule_id)
enrichment_intersect_up <- enrichment(biomolecule_intersect_up, class_and_GO_bp, m3$biomolecule_id)

enrichment_intersect_up[5,]

enrichment_intersect_down_std <- enrichment(EN_merge$standardized_name.COVID[EN_merge$biomolecule_id %in% biomolecule_intersect_down], class_and_GO_bp_standardized_names, EN_merge$standardized_name.COVID)
enrichment_intersect_up_std <- enrichment(EN_merge$standardized_name.COVID[EN_merge$biomolecule_id %in% biomolecule_intersect_up], class_and_GO_bp_standardized_names, EN_merge$standardized_name.COVID)


pdf("plots/Fig2_KAO_GO_BP_Enrichment_based_on_UP_intersect.pdf", width = 7, height = 9)
par(mar = c(4,25,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
par(mfrow = c(2,1))
barplot(-log10(enrichment_intersect_up [order(enrichment_intersect_up[,2]),2][1:15][15:1]), 
        xlab ="-log10(q-value)", 
        cex.names = 0.9, horiz = T, col = 1, 
        names = enrichment_intersect_up [order(enrichment_intersect_up [,2]),1][1:15][15:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significantly up with COVID status and severity")

barplot(-log10(enrichment_intersect_down[order(enrichment_intersect_down[,2]),2][1:15][15:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 3, 
        names = enrichment_intersect_down[order(enrichment_intersect_down[,2]),1][1:15][15:1],
        main = "Enrichment of molecule class and GO Biological processes\nbiomolecules significantly down with COVID status and severity")

dev.off()

write.csv(EN_merge[EN_merge$biomolecule_id %in% biomolecule_interest_directional,], "intersect_features.csv")
# 
# df_forScott<- t_df_wide_exprs[row.names(t_df_wide_exprs) %in% setdiff(biomolecule_interest_directional, biomolecule_unkowns), ]
# colnames(df_forScott) <- df_wide_all$sample_id
# write.csv(df_forScott, "data/Fig2d_feature_table.csv")

write.csv(cbind(enrichment_intersect_up_std, enrichment_intersect_down_std), "data/enrichment_intersect_COVID_and_sevity.csv")
