######## 08_KAO_crossome_correltions.R #######

## Corrlation is very large is all-ome, but if simplified 
## to proteins x metabolite/lipids it might be easier to 
## digest and reveal connections across omes.

library(DBI)
library(RSQLite)
library(pheatmap)
colors <- read.csv("reference/color_palette.txt",  stringsAsFactors = F)[,2]
palette(colors)
source("eda/KAO/0_pathway_toolkit.R")

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
pvalues_HFD <- dbGetQuery(con, "SELECT * FROM pvalues WHERE comparison = 'Hospital_free_days_45'")


# pull
metadata_GO_bp <- dbGetQuery(con, "SELECT * FROM metadata WHERE metadata_type = 'GO_biological_process'")
metadata_GO_mf <- dbGetQuery(con, "SELECT * FROM metadata WHERE metadata_type = 'GO_molecular_function'")
metadata_class <- dbGetQuery(con, "SELECT * FROM metadata WHERE metadata_type = 'Larger_class' OR metadata_type = 'Lipid Class'")

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
protein_hclust <- hclust(dist(cor_4omes_kendall$cor[proteins,metabolites_lipids][filter_row, filter_col]), method = "complete")
plot(protein_hclust, labels = F)
k = 6
protein_hclust_clusters <- cutree(as.hclust(protein_hclust), k=k)

rect.hclust(protein_hclust, k = k, border = "red")

# clustering lipids/metabolites
metabolite_lipid_hclust <- hclust(dist(t(cor_4omes_kendall$cor[proteins,metabolites_lipids][filter_row, filter_col])), method = "complete")
plot(metabolite_lipid_hclust, labels = F)
rect.hclust(metabolite_lipid_hclust, k = 6, border = "red")

metabolite_lipid_hclust_clusters <- cutree(as.hclust(metabolite_lipid_hclust), k = 6)

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
merge_annotation <- merge(pvalues, protein_hclust_clusters,
                          by.x = "biomolecule_id", by.y = 0, all = T)

merge_annotation_2 <- merge(merge_annotation, metabolite_lipid_hclust_clusters, 
                            by.x = "biomolecule_id", by.y = 0, all = T)

merge_annotation_3 <- merge(merge_annotation_2, pvalues_HFD,
                            by.x = "biomolecule_id", by.y = "biomolecule_id", all = T)

names(merge_annotation_3)

annotation_row <- data.frame(protein_clusters = as.factor(merge_annotation_3$y.x), 
                             metabolite_lipid_clusters = as.factor(merge_annotation_3$y.y) ,
                             sig_with_COVID = as.factor(merge_annotation_3$q_value.x < 0.05), 
                             sig_with_HFD = as.factor(merge_annotation_3$q_value.y < 0.05))
row.names(annotation_row) <- merge_annotation_2$biomolecule_id

# annotation colors 
annotation_colors <- list(protein_clusters = colors[c(5:10)], 
                          metabolite_lipid_clusters = colors[c(5:10)], 
                          sig_with_COVID = c("white", colors[1],"white"),
                          sig_with_HFD = c("white", colors[3]))

names(annotation_colors[["protein_clusters"]]) <- as.character(levels(annotation_row$protein_clusters))
names(annotation_colors[["metabolite_lipid_clusters"]]) <- as.character(levels(annotation_row$metabolite_lipid_clusters))
names(annotation_colors[["sig_with_COVID"]]) <- as.character(levels(annotation_row$sig_with_COVID))
names(annotation_colors[["sig_with_HFD"]]) <- as.character(levels(annotation_row$sig_with_HFD))

# Creating heatmap for cross-ome correlation 
pdf("plots/heatmap_cross_ome_correlations_kendall_KAO_v2.pdf", width = 40, height = 30)
pheatmap(cor_4omes_kendall$cor[proteins,metabolites_lipids][filter_row, filter_col],
         annotation_row = annotation_row[c(1,3,4)],
         annotation_col = annotation_row[c(2,3,4)],
         annotation_colors = annotation_colors,
         display_numbers = highlight_sig[filter_row,filter_col],
         labels_col = col_labels,
         labels_row = row_labels,
         cellwidth = 10, cellheight = 10)
dev.off()


# 
# ##### Correlation heatmap: proteins in Row and metabolties in column p < 0.05 #### 
# # Creating a filter for proteins where it must have at least 0.4 Tau with a metabolite
# filter_row <- rowSums(cor_4omes_kendall$adjusted_pvalue[proteins,metabolites_lipids] < 0.05 )>1
# table(filter_row) #152
# 
# # Creating a filter for metabolites-lipids where they must have at least one Tau value over 0.4 
# filter_col <- colSums(cor_4omes_kendall$adjusted_pvalue[proteins,metabolites_lipids] < 0.05 )>1
# table(filter_col) #220
# 
# # row labels as genes 
# row_labels <- df_proteins$geneNames[match(names_cor[proteins][filter_row], df_proteins$biomolecule_id) ]
# 
# df_proteins$biomolecule_id[match(names_cor[proteins][filter_row], df_proteins$biomolecule_id) ]== names_cor[proteins][filter_row]
# 
# # col labels as metabolites ID
# col_labels <- biomolecules$standardized_name[match(names_cor[metabolites_lipids][filter_col], biomolecules$biomolecule_id)]
# 
# #checkmatch
# biomolecules$biomolecule_id[match(names_cor[metabolites_lipids][filter_col], biomolecules$biomolecule_id)] == names_cor[metabolites_lipids][filter_col] 
# 
# # Creating heatmap for cross-ome correlation 
# pdf("heatmap_cross_ome_correlations_kendall_KAO_v3.pdf", width = 50, height = 50)
# pheatmap(cor_4omes_kendall$cor[proteins,metabolites_lipids][filter_row, filter_col],
#          annotation_row = annotation_row,
#          annotation_col = annotation_row,
#          annotation_colors = annotation_colors,
#          display_numbers = highlight_sig[filter_row,filter_col],
#          labels_col = col_labels,
#          labels_row = row_labels,
#          cellwidth = 10, cellheight = 10)
# dev.off()


#### Enrichment analysis and members of interesting clusters ###### 

##### make reference set with go term bp and class info #### 

GO_bp <- make_reference_sets(metadata_GO_bp$metadata_value, metadata_GO_bp$biomolecule_id)
class <- make_reference_sets(metadata_class$metadata_value, metadata_class$biomolecule_id)

class <- class[lapply(class, length) != 1]
GO_bp <- GO_bp[lapply(GO_bp, length) != 1]

### Protein cluster enrichment ##### 
Protein_cluster_1 <- enrichment(names(protein_hclust_clusters)[protein_hclust_clusters == 1], GO_bp, names(protein_hclust_clusters) )
Protein_cluster_2 <- enrichment(names(protein_hclust_clusters)[protein_hclust_clusters == 2], GO_bp, names(protein_hclust_clusters) )
Protein_cluster_3 <- enrichment(names(protein_hclust_clusters)[protein_hclust_clusters == 3], GO_bp, names(protein_hclust_clusters) )
Protein_cluster_4 <- enrichment(names(protein_hclust_clusters)[protein_hclust_clusters == 4], GO_bp, names(protein_hclust_clusters) )
Protein_cluster_5 <- enrichment(names(protein_hclust_clusters)[protein_hclust_clusters == 5], GO_bp, names(protein_hclust_clusters) )
Protein_cluster_6 <- enrichment(names(protein_hclust_clusters)[protein_hclust_clusters == 6], GO_bp, names(protein_hclust_clusters) )

#### Metabolite cluster enrichment #### 

metaboite_cluster_1 <- enrichment(names(metabolite_lipid_hclust_clusters)[metabolite_lipid_hclust_clusters == 1], class, names(metabolite_lipid_hclust_clusters))
metaboite_cluster_2 <- enrichment(names(metabolite_lipid_hclust_clusters)[metabolite_lipid_hclust_clusters == 2], class, names(metabolite_lipid_hclust_clusters))
metaboite_cluster_3 <- enrichment(names(metabolite_lipid_hclust_clusters)[metabolite_lipid_hclust_clusters == 3], class, names(metabolite_lipid_hclust_clusters))
metaboite_cluster_4 <- enrichment(names(metabolite_lipid_hclust_clusters)[metabolite_lipid_hclust_clusters == 4], class, names(metabolite_lipid_hclust_clusters))
metaboite_cluster_5 <- enrichment(names(metabolite_lipid_hclust_clusters)[metabolite_lipid_hclust_clusters == 5], class, names(metabolite_lipid_hclust_clusters))
metaboite_cluster_6 <- enrichment(names(metabolite_lipid_hclust_clusters)[metabolite_lipid_hclust_clusters == 6], class, names(metabolite_lipid_hclust_clusters))

## plot

pdf("plots/08_KAO_enrichment_protein_cluster1_to_4.pdf")
par(mar = c(4,5,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
par(mfrow=c(3,2))
barplot(-log10(Protein_cluster_1[order(Protein_cluster_1[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 5, 
        names = Protein_cluster_1[order(Protein_cluster_1[,2]),1][1:5][5:1],
        main = "Protein_cluster_1")
barplot(-log10(Protein_cluster_2[order(Protein_cluster_2[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 6, 
        names = Protein_cluster_2[order(Protein_cluster_2[,2]),1][1:5][5:1],
        main = "Protein_cluster_2")
barplot(-log10(Protein_cluster_3[order(Protein_cluster_3[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 7, 
        names = Protein_cluster_3[order(Protein_cluster_3[,2]),1][1:5][5:1],
        main = "Protein_cluster_3")
barplot(-log10(Protein_cluster_4[order(Protein_cluster_4[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 8, 
        names = Protein_cluster_4[order(Protein_cluster_4[,2]),1][1:5][5:1],
        main = "Protein_cluster_4")

barplot(-log10(Protein_cluster_5[order(Protein_cluster_5[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 9, 
        names = Protein_cluster_5[order(Protein_cluster_5[,2]),1][1:5][5:1],
        main = "Protein_cluster_5")

barplot(-log10(Protein_cluster_6[order(Protein_cluster_6[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 10, 
        names = Protein_cluster_6[order(Protein_cluster_6[,2]),1][1:5][5:1],
        main = "Protein_cluster_6")
dev.off()

#### ploting metabolite clusters ####

pdf("plots/08_KAO_enrichment_metabolite_cluster1_to_4.pdf")
par(mar = c(4,5,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
par(mfrow=c(3,2))
barplot(-log10(metaboite_cluster_1[order(metaboite_cluster_1[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 5, 
        names = metaboite_cluster_1[order(metaboite_cluster_1[,2]),1][1:5][5:1],
        main = "metabolite_cluster_1")
barplot(-log10(metaboite_cluster_2[order(metaboite_cluster_2[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 6, 
        names = metaboite_cluster_2[order(metaboite_cluster_2[,2]),1][1:5][5:1],
        main = "metaboite_cluster_2")
barplot(-log10(metaboite_cluster_3[order(metaboite_cluster_3[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 7, 
        names = metaboite_cluster_3[order(metaboite_cluster_3[,2]),1][1:5][5:1],
        main = "metaboite_cluster_3")
barplot(-log10(metaboite_cluster_4[order(metaboite_cluster_4[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 8, 
        names = metaboite_cluster_4[order(metaboite_cluster_4[,2]),1][1:5][5:1],
        main = "metaboite_cluster_4")

barplot(-log10(metaboite_cluster_5[order(metaboite_cluster_5[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 9, 
        names = metaboite_cluster_5[order(metaboite_cluster_5[,2]),1][1:5][5:1],
        main = "metaboite_cluster_5")

barplot(-log10(metaboite_cluster_6[order(metaboite_cluster_6[,2]),2][1:5][5:1]), 
        xlab ="-log10(p-value)", 
        cex.names = 0.9, horiz = T, col = 10, 
        names = metaboite_cluster_6[order(metaboite_cluster_6[,2]),1][1:5][5:1],
        main = "metaboite_cluster_6")
dev.off()



# Creating heatmap for cross-ome correlation 
pdf("plots/heatmap_cross_ome_correlations_kendall_KAO_v3.pdf", width = 40, height = 30)
pheatmap(cor_4omes_kendall$cor[proteins,metabolites_lipids][filter_row, filter_col],
         annotation_row = annotation_row[c(1,3,4)],
         annotation_col = annotation_row[c(2,3,4)],
         annotation_colors = annotation_colors,
         display_numbers = highlight_sig[filter_row,filter_col],
         show_colnames = F, show_rownames = F,
         #labels_col = col_labels,
         #labels_row = row_labels,
         cellwidth = 2, cellheight = 2)
dev.off()


#### Zoom in heatmap of interesting clusters #### 

protein_clusters_of_interest <- protein_hclust_clusters == 2 | protein_hclust_clusters == 5 | protein_hclust_clusters == 6
metabolite_clusters_of_interest <- metabolite_lipid_hclust_clusters == 1 | metabolite_lipid_hclust_clusters == 2

pdf("plots/heatmap_cross_ome_correlations_kendall_KAO_v4.pdf", width = 15, height = 15)
pheatmap(cor_4omes_kendall$cor[proteins,metabolites_lipids][filter_row, filter_col][protein_clusters_of_interest, metabolite_clusters_of_interest],
         annotation_row = annotation_row[c(1,3,4)],
         annotation_col = annotation_row[c(2,3,4)],
         annotation_colors = annotation_colors,
         display_numbers = highlight_sig[filter_row,filter_col][protein_clusters_of_interest, metabolite_clusters_of_interest],
         labels_col = col_labels[metabolite_clusters_of_interest],
         labels_row = row_labels[protein_clusters_of_interest],
         cellwidth = 10, cellheight = 10)
dev.off()


