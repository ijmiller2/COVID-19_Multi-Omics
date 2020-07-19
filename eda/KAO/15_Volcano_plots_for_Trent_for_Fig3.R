##### 15_Volcano_plots_for_Trent_for_Fig3.R ##### 

## These data pull infromation from pvalues table and provides numbers for the text about figure 2. 

library(DBI)
library(RSQLite)

source("eda/KAO/0_pathway_toolkit.R")
library(pheatmap)

colors <- read.csv("reference/color_palette.txt",  stringsAsFactors = F)[,2]
palette(colors)

### Fig 2A - volcano plot COVID. #### 

# connect 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# pull pvalues 
pvalues_COVID_proteins <- dbGetQuery(con, "SELECT * FROM pvalues 
          INNER JOIN biomolecules on biomolecules.biomolecule_id = pvalues.biomolecule_id
          INNER JOIN metadata on metadata.biomolecule_id = pvalues.biomolecule_id
          WHERE pvalues.comparison = 'COVID_vs_NONCOVID' AND biomolecules.keep = 1 
          AND biomolecules.omics_id = 1 AND metadata.metadata_type = 'gene_name' 
                  ")

# disconnect
dbDisconnect(con) 

## load in log2 fold change values 

FC <- read.csv("data/COVID_fc_by_biomolecule_ID.csv")

## load in Trent's protein groups 

fig3 <- read.csv("data/Proteins grouped for Fig 3 Volcano Plots.csv", stringsAsFactors = F)
fig3 <- fig3[,1:2]
## merge pvalues and FC by biomolecule_id

pvalues_COVID_proteins <- pvalues_COVID_proteins[,-grep('biomolecule_id', names(pvalues_COVID_proteins))[2:3]]

m <- merge(pvalues_COVID_proteins, FC, by.x = 'biomolecule_id', by.y = 'biomolecule_id')
m <- data.frame(m$biomolecule_id, m$q_value, m$FC, m$metadata_value)

geneName <- apply(m, 1, function(x) strsplit(x[4], ';')[[1]][1])
m$geneName <- geneName

fig3$X.1[122] <- "APOB"
fig3$X.1[30] <- "ITGA2B"
m2 <- merge(m, fig3, by.x = 'geneName', by.y = 'X.1', all = T)

unique(m2$X)

pdf("plots/ForTrent_Fig3_KAO_20200626.pdf", useDingbats = F)
par(mfrow = c(3,2))

for(i in unique(m2$X)[-1]) {

  plot(m2$m.FC, -log(m2$m.q_value), 
       pch = 19, col = "gray70",
       ylim = c(0, 30), xlim = c(-3, 4),
       main = paste(i), 
       ylab = "-log(adjusted p-value of likelyhood ratio)",
       xlab = "log2(abundance COVID-19 positive/COVID-19 negative)",
       las = 1, bty = "l", cex.axis = 1.2)
  points(m2$m.FC[m2$X == i], -log(m2$m.q_value)[m2$X == i], 
         pch = 19, col = "gray10")
 
}
dev.off()



## These data pull infromation from pvalues table and provides numbers for the text about figure 2. 

library(DBI)
library(RSQLite)

source("eda/KAO/0_pathway_toolkit.R")
library(pheatmap)
library(ellipse)
library(scales)
colors <- read.csv("reference/color_palette.txt",  stringsAsFactors = F)[,2]
palette(colors)

### Fig 3 - volcano plot COVID-status #### 

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


pdf("plots/volcano_plots_all_omes_fig3.pdf", useDingbats = F, width = 9, height = 1.5)

par(mfrow = c(1,6), mar = c(3,4,4,1), cex = 0.4)
help(par)
for(i in c(118)){
  up <- class_and_GO_bp[[enrichment_significant_up_COVID[i,1]]]
  
  plot(m$FC[!unknown], -log(m$q_value)[!unknown], 
       pch = 19, col = "gray70",
       ylim = c(0, 25), xlim = c(-4, 6),
       main = paste(enrichment_significant_up_COVID[i,1]), 
       ylab = "-log(adjusted p-value of likelyhood ratio)",
       xlab = "log2(abundance COVID-19 positive/COVID-19 negative)",
       las = 1, bty = "l", cex.axis = 1.2)
  points(m$FC[m$biomolecule_id %in% up], -log(m$q_value)[m$biomolecule_id %in% up], 
         pch = 19, col = colors[c(12,9,10,11,11)][m$omics_id[m$biomolecule_id %in% up]])
  
}


for(i in c(313)){
  up <- class_and_GO_bp[[enrichment_significant_up_COVID[i,1]]]
  
  plot(m$FC[!unknown], -log(m$q_value)[!unknown], 
       pch = 19, col = "gray70",
       ylim = c(0, 25), xlim = c(-4, 6),
       main = paste(enrichment_significant_up_COVID[i,1]), 
       ylab = "-log(adjusted p-value of likelyhood ratio)",
       xlab = "log2(abundance COVID-19 positive/COVID-19 negative)",
       las = 1, bty = "l", cex.axis = 1.2)
  points(m$FC[m$biomolecule_id %in% up], -log(m$q_value)[m$biomolecule_id %in% up], 
         pch = 19, col = colors[c(12,9,10,11,11)][m$omics_id[m$biomolecule_id %in% up]])
  
}


for(i in c(4)){
  up <- class_and_GO_bp[[enrichment_significant_up_COVID[i,1]]]
  
  plot(m$FC[!unknown], -log(m$q_value)[!unknown], 
       pch = 19, col = "gray70",
       ylim = c(0, 25), xlim = c(-4, 6),
       main = paste(enrichment_significant_up_COVID[i,1]), 
       ylab = "-log(adjusted p-value of likelyhood ratio)",
       xlab = "log2(abundance COVID-19 positive/COVID-19 negative)",
       las = 1, bty = "l", cex.axis = 1.2)
  points(m$FC[m$biomolecule_id %in% up], -log(m$q_value)[m$biomolecule_id %in% up], 
         pch = 19, col = colors[c(12,9,10,11,11)][m$omics_id[m$biomolecule_id %in% up]])
  
}

for(i in c(28)){
  up <- class_and_GO_bp[[enrichment_significant_up_COVID[i,1]]]
  
  plot(m$FC[!unknown], -log(m$q_value)[!unknown], 
       pch = 19, col = "gray70",
       ylim = c(0, 25), xlim = c(-4, 6),
       main = paste(enrichment_significant_up_COVID[i,1]), 
       ylab = "-log(adjusted p-value of likelyhood ratio)",
       xlab = "log2(abundance COVID-19 positive/COVID-19 negative)",
       las = 1, bty = "l", cex.axis = 1.2)
  points(m$FC[m$biomolecule_id %in% up], -log(m$q_value)[m$biomolecule_id %in% up], 
         pch = 19, col = colors[c(12,9,10,11,11)][m$omics_id[m$biomolecule_id %in% up]])
  
}


for(i in c(56)){
  up <- class_and_GO_bp[[enrichment_significant_up_COVID[i,1]]]
  
  plot(m$FC[!unknown], -log(m$q_value)[!unknown], 
       pch = 19, col = "gray70",
       ylim = c(0, 25), xlim = c(-4, 6),
       main = paste(enrichment_significant_up_COVID[i,1]), 
       ylab = "-log(adjusted p-value of likelyhood ratio)",
       xlab = "log2(abundance COVID-19 positive/COVID-19 negative)",
       las = 1, bty = "l", cex.axis = 1.2)
  points(m$FC[m$biomolecule_id %in% up], -log(m$q_value)[m$biomolecule_id %in% up], 
         pch = 19, col = colors[c(12,9,10,11,11)][m$omics_id[m$biomolecule_id %in% up]])
  
}


for(i in c(132)){
  up <- class_and_GO_bp[[enrichment_significant_up_COVID[i,1]]]
  
  plot(m$FC[!unknown], -log(m$q_value)[!unknown], 
       pch = 19, col = "gray70",
       ylim = c(0, 35), xlim = c(-4, 6),
       main = paste(enrichment_significant_up_COVID[i,1]), 
       ylab = "-log(adjusted p-value of likelyhood ratio)",
       xlab = "log2(abundance COVID-19 positive/COVID-19 negative)",
       las = 1, bty = "l", cex.axis = 1.2)
  points(m$FC[m$biomolecule_id %in% up], -log(m$q_value)[m$biomolecule_id %in% up], 
         pch = 19,col = colors[c(12,9,10,11,11)][m$omics_id[m$biomolecule_id %in% up]])
  
}

dev.off()

## source("14_KAO_figure_2_version_1.R")
cex_breaks = seq(0,35, length.out = 5)
cex_bin_values <- .bincode(-log(EN_merge$q_value.HFD), cex_breaks, include.lowest = T)
cex_HFD <- seq(0.5, 3, length.out = 5)[cex_bin_values]
cex_sig <- c(0.7, 1.2)[as.numeric(EN_merge$biomolecule_id %in% biomolecule_interest_directional)+1]

unknown <- grepl("nknown", EN_merge$standardized_name.COVID)

pdf("plots/Volcano_plot_neutrophil_degranulation_v2.pdf", height = 5, width = 5, useDingbats = F)
for(i in c(132)){
  up <- class_and_GO_bp[[enrichment_significant_up_COVID[i,1]]]
  par(mgp = c(1.8, 0.7, 0), tcl = -0.3)
  
  plot(EN_merge$FC[!unknown], -log(EN_merge$q_value.COVID)[!unknown], 
       pch = 19, col = "gray70",
       cex = 0.7,
       ylim = c(0, 35), xlim = c(-4, 6),
       main = paste(enrichment_significant_up_COVID[i,1]), 
       ylab = "-log(FDR)",
       xlab = "log2(abundance COVID-19 positive/COVID-19 negative)",
       las = 1, bty = "l", cex.axis = 1.2)
  points(EN_merge$FC[EN_merge$biomolecule_id %in% up], -log(EN_merge$q_value.COVID)[EN_merge$biomolecule_id %in% up], 
         pch = 19,
         cex = cex_sig[EN_merge$biomolecule_id %in% up],
         col = colors[c(12,9,10,11,11)][EN_merge$omics_id[EN_merge$biomolecule_id %in% up]])
  text(EN_merge$FC[(EN_merge$biomolecule_id %in% up) & (EN_merge$biomolecule_id %in% biomolecule_interest_directional)], 
       -log(EN_merge$q_value.COVID)[(EN_merge$biomolecule_id %in% up) & (EN_merge$biomolecule_id %in% biomolecule_interest_directional)],
       EN_merge$standardized_name.COVID[(EN_merge$biomolecule_id %in% up) & (EN_merge$biomolecule_id %in% biomolecule_interest_directional)],
       pos =4
  )

}
dev.off()

enrichment_significant_up_COVID_protein <- enrichment(m$standardized_name[m$q_value < 0.05 & m$FC > 0 & m$omics_id == 1], class_and_GO_bp_standardized_names, m$standardized_name)
enrichment_significant_down_COVID_protein <- enrichment(m$standardized_name[m$q_value < 0.05 & m$FC < 0 & m$omics_id == 1], class_and_GO_bp_standardized_names, m$standardized_name)

enrichment_significant_up_COVID_transcipts <- enrichment(m$standardized_name[m$q_value < 0.05 & m$FC > 0 & m$omics_id == 5], class_and_GO_bp_standardized_names, m$standardized_name)
enrichment_significant_down_COVID_transcripts <- enrichment(m$standardized_name[m$q_value < 0.05 & m$FC < 0 & m$omics_id == 5], class_and_GO_bp_standardized_names, m$standardized_name)

write.csv(cbind(enrichment_significant_up_COVID_protein, enrichment_significant_down_COVID_protein, enrichment_significant_up_COVID_transcipts, enrichment_significant_down_COVID_transcripts), file = "data/differences_COVID_enrichment_proteins_transcipts.csv")
