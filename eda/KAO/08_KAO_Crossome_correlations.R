######## 08_KAO_crossome_correltions.R #######

## Corrlation is very large is all-ome, but if simplified 
## to proteins x metabolite/lipids it might be easier to 
## digest and reveal connections across omes.

library(DBI)
library(RSQLite)
library(pheatmap)
colors <- read.csv("reference/color_palette.txt",  stringsAsFactors = F)[,2]
palette(colors)

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

pvalues <- dbGetQuery(con, "SELECT * FROM pvalues WHERE comparison = 'COVID_vs_NONCOVID'")

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

# Creating a filter for proteins where it must have at least 0.4 Tau with a metabolite
filter_row <- rowSums(abs(cor_4omes_kendall$cor[proteins,metabolites_lipids]) > 0.4 )>1
table(filter_row) #112

# Creating a filter for metabolites-lipids where they must have at least one Tau value over 0.4 
filter_col <- colSums(abs(cor_4omes_kendall$cor[proteins,metabolites_lipids]) > 0.4 )>1
table(filter_col) #118

# significance star
highlight_sig <- cor_4omes_kendall[[3]][proteins,metabolites_lipids] <  0.05
highlight_sig[!highlight_sig] <- ""
highlight_sig[highlight_sig == TRUE] <- "*"


###### Extracting out cluster information #########
# clustering proteins
protein_hclust <- hclust(dist(t(cor_4omes_kendall$cor[proteins,metabolites_lipids][filter_row, filter_col])), method = "complete")
plot(protein_hclust, labels = F)
k = 8
protein_hclust_clusters <- cutree(as.hclust(protein_hclust), k=k)

rect.hclust(protein_hclust, k = k, border = "red")

clusplot(t(cor_4omes_kendall$cor[proteins,metabolites_lipids][filter_row, filter_col]), protein_hclust_clusters, color =T, shade = T, lines = 0, col.clus = 1:10)

# clustering lipids/metabolites
metabolite_lipid_hclust <- hclust(dist(cor_4omes_kendall$cor[proteins,metabolites_lipids][filter_row, filter_col]), method = "complete")
plot(metabolite_lipid_hclust, labels = F)
rect.hclust(metabolite_lipid_hclust, k = 7, border = "red")

metabolite_lipid_hclust_clusters <- cutree(as.hclust(metabolite_lipid_hclust), k = 7)

###### Annotations for heatmap #######

# for plotting,extract gene names from proteins-metadata table.
df_proteins$geneNames<- apply(df_proteins , 1, function(x) strsplit(strsplit(x[4], "GN=")[[1]][2], " ")[[1]][1])

# row labels as genes 
row_labels <- df_proteins$geneNames[match(names_cor[proteins][filter_row], df_proteins$biomolecule_id) ]

df_proteins$biomolecule_id[match(names_cor[proteins][filter_row], df_proteins$biomolecule_id) ]== names_cor[proteins][filter_row]

# col labels as metabolites ID
col_labels <- biomolecules$standardized_name[match(names_cor[metabolites_lipids][filter_col], biomolecules$biomolecule_id)]

#checkmatch
biomolecules$biomolecule_id[match(names_cor[metabolites_lipids][filter_col], biomolecules$biomolecule_id)] == names_cor[metabolites_lipids][filter_col] 

# annotation based on significant with COVID and cluster member
merge_annotation <- merge(pvalues, protein_hclust_clusters, by.x = "biomolecule_id", by.y = 0)
merge_annotation_2 <- merge(merge_annotation, metabolite_lipid_hclust_clusters, by.x = "biomolecule.id", by.y = 0)

annotation_row <- data.frame(protein_clusters = as.factor(merge_annotation_row$protein_hclust_clusters), metabolite_lipid_clusters = as.factor(merge_annotation_2$metabolite_lipid_hclust_clusters) ,sig_with_COVID = as.factor(merge_annotation_row$q_value < 0.05))
row.names(annotation_row) <- merge_annotation_row$biomolecule_id

# annotation colors 
annotation_colors <- list(protein_clusters = c(5:12), metabolite_lipid_clusters = c(5:11), sig_with_COVID = c(3,1,"white"))

names(annotation_colors[["protein_clusters"]]) <- as.character(levels(annotation_row$protein_clusters))
names(annotation_colors[["metabolite_lipid_clusters"]]) <- as.character(levels(annotation_row$metabolite_lipid_clusters))
names(annotation_colors[["sig_with_COVID"]]) <- as.character(levels(annotation_row$sig_with_COVID))

# Creating heatmap for cross-ome correlation 
pdf("heatmap_cross_ome_correlations_kendall_KAO_v2.pdf", width = 40, height = 30)
pheatmap(cor_4omes_kendall$cor[proteins,metabolites_lipids][filter_row, filter_col],
         annotation_row = annotation_row[c(1,3)],
         annotation_col = annotation_row[c(2,3)],
         annotation_colors = annotation_colors,
         display_numbers = highlight_sig[filter_row,filter_col],
         labels_col = col_labels,
         labels_row = row_labels,
         cellwidth = 10, cellheight = 10)
dev.off()

##### Corrlation heatmap: proteins in Row and metabolties in column p < 0.05 #### 
# Creating a filter for proteins where it must have at least 0.4 Tau with a metabolite
filter_row <- rowSums(cor_4omes_kendall$adjusted_pvalue[proteins,metabolites_lipids] < 0.05 )>1
table(filter_row) #152

# Creating a filter for metabolites-lipids where they must have at least one Tau value over 0.4 
filter_col <- colSums(cor_4omes_kendall$adjusted_pvalue[proteins,metabolites_lipids] < 0.05 )>1
table(filter_col) #220

# row labels as genes 
row_labels <- df_proteins$geneNames[match(names_cor[proteins][filter_row], df_proteins$biomolecule_id) ]

df_proteins$biomolecule_id[match(names_cor[proteins][filter_row], df_proteins$biomolecule_id) ]== names_cor[proteins][filter_row]

# col labels as metabolites ID
col_labels <- biomolecules$standardized_name[match(names_cor[metabolites_lipids][filter_col], biomolecules$biomolecule_id)]

#checkmatch
biomolecules$biomolecule_id[match(names_cor[metabolites_lipids][filter_col], biomolecules$biomolecule_id)] == names_cor[metabolites_lipids][filter_col] 

# Creating heatmap for cross-ome correlation 
pdf("heatmap_cross_ome_correlations_kendall_KAO_v3.pdf", width = 50, height = 50)
pheatmap(cor_4omes_kendall$cor[proteins,metabolites_lipids][filter_row, filter_col],
         annotation_row = annotation_row,
         annotation_col = annotation_row,
         annotation_colors = annotation_colors,
         display_numbers = highlight_sig[filter_row,filter_col],
         labels_col = col_labels,
         labels_row = row_labels,
         cellwidth = 10, cellheight = 10)
dev.off()
