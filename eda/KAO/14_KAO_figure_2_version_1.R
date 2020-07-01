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
metadata_GO_bp <- dbGetQuery(con, "SELECT * FROM metadata WHERE metadata_type = 'GO_biological_process'")
metadata_class <- dbGetQuery(con, "SELECT * FROM metadata WHERE metadata_type = 'Larger_class' OR metadata_type = 'Lipid Class'")

# disconnect
dbDisconnect(con) 

# 

metadata_class <- metadata_class[!(metadata_class$metadata_value == ""),]
metadata_class$metadata_value <- paste(metadata_class$metadata_value,"[Metabolite Class]")


##### make reference set with go term bp and class info #### 
metadata_GO_bp_class<- rbind(metadata_GO_bp, metadata_class)

class_and_GO_bp <- make_reference_sets(metadata_GO_bp_class$metadata_value, metadata_GO_bp_class$biomolecule_id)

#### performing enrichment based LR test significant for COVID #### 
enrichment_significant_up_COVID <- enrichment(m$biomolecule_id[m$q_value < 0.05 & m$FC > 1], class_and_GO_bp, m$biomolecule_id)
enrichment_significant_down_COVID <- enrichment(m$biomolecule_id[m$q_value < 0.05 & m$FC < -1], class_and_GO_bp, m$biomolecule_id)


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

both <- m2$biomolecule_id[m2$q_value.x < 0.05 & m2$q_value.y < 0.05] 

names(df_wide_exprs) <- sub("normalized_abundance.", "", names(df_wide_exprs)) 
t_df_wide_exprs <- t(df_wide_exprs)

table( row.names(t_df_wide_exprs) %in% both)

m2$standardized_name.y[m2$q_value.x < 0.05 & m2$q_value.y < 0.05] 
## heatmap 

pheatmap(t_df_wide_exprs[row.names(t_df_wide_exprs) %in% both, ], scale = "row")
