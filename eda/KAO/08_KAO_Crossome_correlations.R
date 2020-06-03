######## 08_KAO_crossome_correltions.R #######

## Corrlation is very large is all-ome, but if simplified 
## to proteins x metabolite/lipids it might be easier to 
## digest and reveal connections across omes.

library(DBI)
library(RSQLite)
library(pheatmap)

###### Loading in data from Yuchen, kendall correlations #####

load("P:/All_20200428_COVID_plasma_multiomics/Correlation/cor_4omes_kendall.RData")

#### Establish a connection to the DB #####
con <- dbConnect(RSQLite::SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")

#### Pull data from DB ####

dbListTables(con)

df_proteins<- dbGetQuery(con, "SELECT * FROM metadata
            WHERE metadata_type = 'fasta_header'  
           ")

biomolecules <-dbGetQuery(con, "SELECT * FROM biomolecules")
dbDisconnect(con)

##### Find gelsolin biomolecule_id ####### 
gelsolin <- grep("Gelsolin", df_proteins$metadata_value)

df_proteins[gelsolin, ] # biomolecule_id = 7974

##### Corrlation heatmap: proteins in Row and metabolties in column#### 
#names holds biomolecule_id
names_cor <- row.names(cor_4omes_kendall$cor)

# matching biomolecule to determine ome, also filtering keep =1 for updated lipids filter
proteins <- names_cor %in% biomolecules$biomolecule_id[biomolecules$omics_id==1]
metabolites_lipids <- names_cor %in% biomolecules$biomolecule_id[biomolecules$keep == 1] & !proteins

table(proteins) #517 proteins
table(metabolites_lipids) #3512 metabolites and lipids 

# Creating a filter for proteins where it must have at least one Tau value over 0.4 with a metabolite
filter_row <- rowSums(abs(cor_4omes_kendall$cor[proteins,metabolites_lipids]) > 0.4)>1
table(filter_row) #112 

# Creating a filter for metabolites-lipids where they must have at least one Tau value over 0.4 
filter_col <- colSums(abs(cor_4omes_kendall$cor[proteins,metabolites_lipids]) > 0.4)>1
table(filter_col) #118

# for plotting,extract gene names from proteins-metadata table.
geneNames<- apply(df_proteins , 1, function(x) strsplit(strsplit(x[4], "GN=")[[1]][2], " ")[[1]][1])
  row_labels <- geneNames[df_proteins$biomolecule_id %in% row.names(cor_4omes_kendall$cor)[proteins]][filter_row]

# Creating heatmap for cross-ome correlation 
pdf("heatmap_cross_ome_correlations_kendall_KAO.pdf", width = 40, height = 40)
pheatmap(cor_4omes_kendall$cor[proteins,metabolites_lipids][filter_row, filter_col],
         labels_col = biomolecules$standardized_name[biomolecules$biomolecule_id %in% row.names(cor_4omes_kendall$cor)[metabolites_lipids][filter_col]],
         labels_row = row_labels, cellwidth = 10, cellheight = 10)
dev.off()
