##### 14_KAO_figure_2_version_1.R ##### 

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
pvalues_COVID <- dbGetQuery(con, "SELECT * FROM pvalues 
          INNER JOIN biomolecules on biomolecules.biomolecule_id = pvalues.biomolecule_id
          WHERE pvalues.comparison = 'COVID_vs_NONCOVID' AND biomolecules.keep = 1
                  ")

# disconnect
dbDisconnect(con) 

## load in log2 fold change values 

FC <- read.csv("data/COVID_fc_by_biomolecule_ID.csv")

## merge pvalues and FC by biomolecule_id

pvalues_COVID <- pvalues_COVID[,-grep('biomolecule_id', names(pvalues_COVID))[2]]

m <- merge(pvalues_COVID, FC, by.x = 'biomolecule_id', by.y = 'biomolecule_id')
m$omics_id <- as.factor(m$omics_id)
m$omics_id <- relevel(m$omics_id, ref = 5)
m <- m[order(m$omics_id), ]

unknown <- grepl("nknown", m$standardized_name)
table(unknown)

## Plot Volcano plot 
pdf("plots/KAO_fig_2A_Volcano_plot_v1.pdf", useDingbats = F)

plot(m$FC[unknown], -log(m$q_value)[unknown],
     pch = 1, col = "gray40",
     ylim = c(0, 40), xlim = c(-4, 6),
     main = "Effect of COVID", 
     ylab = "-log(adjusted p-value of likelyhood ratio)",
     xlab = "log2(abundance COVID-19 positive/COVID-19 negative)",
     las = 1, bty = "l", cex.axis = 1.2)
points(m$FC[!unknown], -log(m$q_value)[!unknown], 
     pch = 19, col = colors[c(12,9,10,11,11)][m$omics_id[!unknown]])
     
legend("topleft", c("transcripts", "proteins", "lipids", "small molecules", "unidentified features"), 
       col = c(colors[c(12,9,10,11)], "gray40"), pch = c(19,19,19,19,1),
       cex= 0.8, bty = "n")

dev.off()

## for text:

table(m$omics_id[!unknown],m$q_value[!unknown] < 0.05)

# FALSE TRUE
# 5  6764 6499
# 1   371  146
# 2   478  168
# 3    58    8
# 4    34    5

table(m$q_value[unknown] < 0.05)

# FALSE  TRUE 
# 2245   511 

##### Gene set enrichment 

## Establish a connection to the DB 
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

# pull
metadata_GO_bp <- dbGetQuery(con, "SELECT * FROM metadata WHERE metadata_type = 'GO_biological_process'")
metadata_class <- dbGetQuery(con, "SELECT * FROM metadata WHERE metadata_type = 'Larger_class' OR metadata_type = 'Lipid Class'")

# disconnect
dbDisconnect(con) 

##### make reference set with go term bp and class info #### 
metadata_GO_bp_class<- rbind(metadata_GO_bp, metadata_class)

class_and_GO_bp <- make_reference_sets(metadata_GO_bp_class$metadata_value, metadata_GO_bp_class$biomolecule_id)

#### performing enrichment based LR test significant for COVID #### 
enrichment_significant_up_COVID <- enrichment(m$biomolecule_id[m$q_value < 0.05 & m$FC > 1], class_and_GO_bp, m$biomolecule_id)
enrichment_significant_down_COVID <- enrichment(m$biomolecule_id[m$q_value < 0.05 & m$FC < -1], class_and_GO_bp, m$biomolecule_id)

enrichment_significant_up_COVID[order(enrichment_significant_up_COVID$pvalue),][1,1]
enrichment_significant_down_COVID[order(enrichment_significant_down_COVID$pvalue),][1:6,]

pdf("plots/enriched_in_covid_volcano_plots.pdf", useDingbats = )
par(mfrow = c(3,2))
for(i in 1:3){
up <- class_and_GO_bp[[enrichment_significant_up_COVID[order(enrichment_significant_up_COVID$pvalue),][i,1]]]

plot(m$FC[!unknown], -log(m$q_value)[!unknown], 
     pch = 19, col = "gray70",
     ylim = c(0, 40), xlim = c(-4, 6),
     main = paste(enrichment_significant_up_COVID[order(enrichment_significant_up_COVID$pvalue),][i,1]), 
     ylab = "-log(adjusted p-value of likelyhood ratio)",
     xlab = "log2(abundance COVID-19 positive/COVID-19 negative)",
     las = 1, bty = "l", cex.axis = 1.2)
points(m$FC[m$biomolecule_id %in% up], -log(m$q_value)[m$biomolecule_id %in% up], 
       pch = 19, col = "black")

}

for(i in 1:3){
  down <- class_and_GO_bp[[enrichment_significant_down_COVID[order(enrichment_significant_down_COVID$pvalue),][i,1]]]
  
  plot(m$FC[!unknown], -log(m$q_value)[!unknown], 
       pch = 19, col = "gray70",
       ylim = c(0, 40), xlim = c(-4, 6),
       main = paste(enrichment_significant_down_COVID[order(enrichment_significant_down_COVID$pvalue),][i,1]), 
       ylab = "-log(adjusted p-value of likelyhood ratio)",
       xlab = "log2(abundance COVID-19 positive/COVID-19 negative)",
       las = 1, bty = "l", cex.axis = 1.2)
  points(m$FC[m$biomolecule_id %in% down], -log(m$q_value)[m$biomolecule_id %in% down], 
         pch = 19, col = "black")
  
}
dev.off()