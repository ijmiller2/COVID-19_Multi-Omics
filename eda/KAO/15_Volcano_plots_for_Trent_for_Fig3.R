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
