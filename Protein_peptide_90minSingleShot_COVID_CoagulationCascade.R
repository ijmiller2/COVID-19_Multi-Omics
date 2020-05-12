library(pcaMethods)
library(ggplot2)
library(plotly)
library(reshape2)
library(plyr)
library(pheatmap)

#### Load data ####

plasma_proteinGroups_coagulationCascade <- read.csv("P:/All_20200428_COVID_plasma_multiomics/Proteomics/90 min single-shot reps for Anji/combined_includingIsoforms_Deane_coagulationPathway/txt/peptides_Plasma_90SingleShot.csv", 
                                           header = TRUE, sep = ",", stringsAsFactors = FALSE)

plasma_peptide_coagulationCascade <- read.csv("P:/All_20200428_COVID_plasma_multiomics/Proteomics/90 min single-shot reps for Anji/combined_includingIsoforms_Deane_coagulationPathway/txt/peptides_Plasma_90SingleShot.csv", 
                                                    header = TRUE, sep = ",", stringsAsFactors = FALSE)

uniprot_list <- read.csv("P:/All_20200428_COVID_plasma_multiomics/Proteomics/90 min single-shot reps for Anji/UniprotDatabase_CoagulationCascade_List.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)



# subset specific columns and bind function
subsetLFQ <- function(x){
  y<- x[,which(names(x) %in% c("Protein.IDs", 
                               "Majority.protein.IDs",
                               "Gene.names",
                               "Fasta.headers",
                               "Peptides",
                               "Razor...unique.peptides",
                               "Unique.peptides",
                               "Sequence.coverage....",
                               "Unique.sequence.coverage....",
                               "Mol..weight..kDa.",
                               "Q.value",
                               "Score",
                               "Intensity",
                               "MS.MS.count",
                               "id"))]
  z <- x[,grep("LFQ.intensity.",names(x))]
  z[z == 0] <- NA
  x <- cbind(y,z)
  return(x)
}
# subset columns
subset_proteins_CoagCascade <- subsetLFQ(plasma_proteinGroups_coagulationCascade)

# use id as identifier in row names
rownames(subset_proteins_CoagCascade) <- subset_proteins_CoagCascade$id

# Isolate only the first protein in protein group for Majority.protein.IDs and Protein.IDs
majorityProtein <- lapply(strsplit(subset_proteins_CoagCascade$Majority.protein.IDs, ";"), '[', 1)
protein.ID <- lapply(strsplit(subset_proteins_CoagCascade$Protein.IDs, ";"), '[', 1)

# Replace with reduced array
subset_proteins_CoagCascade$Protein.IDs <- protein.ID
subset_proteins_CoagCascade$Majority.protein.IDs <- majorityProtein

#collapse first two columns
subset_proteins_CoagCascade$Protein.IDs <- vapply(subset_proteins_CoagCascade$Protein.IDs, paste, collapse = ", ", character(1L))
subset_proteins_CoagCascade$Majority.protein.IDs <- vapply(subset_proteins_CoagCascade$Majority.protein.IDs, paste, collapse = ", ", character(1L))

#Function for 50% Filter
remove.features.50percentcuttoff <- function (x) {
  
  # Calculate how many missing values per feature
  features.missing = rowMeans(is.na(x)) 
  print(paste0("Number of protein groups that have over 50% missing measurements: ",sum(features.missing > 0.50))) 
  features.missing.50more = rownames(x)[features.missing > 0.50] 
  
  # create a vector of the features that meet criteria
  keep.features = rownames(x)[features.missing <= 0.50]
  print(paste0("Protein groups that pass the 50% filteration: ", length(keep.features)))
  #names(keep.features) = keep.features 
  
  # create a vector of the features that will be removed from data.frame
  remove.features = rownames(x)[features.missing > 0.50]
  print(paste0("Number of protein groups removed from dataset: ", length(remove.features)))
  #names(remove.features) = remove.features
  
  # If there isn't any features to remove keep all features
  if(sum(features.missing) == 0 ){
    filtered = x[which(rownames(x) %in% keep.features),]
    # otherwise filter out the features that contain over 50% of data missing
  } else{
    filtered = x[-which(rownames(x) %in% remove.features),]
  }
  return(filtered)
}

#How many zero values?
res <- colSums(subset_proteins_CoagCascade[,15:17]==0)/nrow(subset_proteins_CoagCascade[,15:17])*100

# Perform feature removal on only the LFQ data only
clean_df <- remove.features.50percentcuttoff(subset_proteins_CoagCascade[16:ncol(subset_proteins_CoagCascade)])
# Remove Intensity from column names
colnames(clean_df) <- sub(".*intensity.", "", colnames(clean_df))

# Log2 clean_df
clean_df_log2 <- log2(clean_df)

# Create a meta object with only the meta data from MaxQuant
# There aren't any zero values in dataset no need to create a filtered meta data
meta <- subset_proteins_CoagCascade[,1:15]

#### Plot Protein Intensity histograms ####
hist(clean_df_log2$plasma_2.2_rep1,
     breaks = 15,
     ylim = c(0,7),
     xlim = c(26,40),
     xlab = "Log2(LFQ)",
     main = "Histogram 90 min Log2(LFQ) COVID19 Sample1")

hist(clean_df_log2$plasma_2.2_rep2,
     breaks = 15,
     ylim = c(0,7),
     xlim = c(26,40),
     xlab = "Log2(LFQ)",
     main = "Histogram 90 min Log2(LFQ) COVID19 Sample2")

hist(clean_df_log2$plasma_2.2_rep3,
     breaks = 15,
     ylim = c(0,7),
     xlim = c(26,40),
     xlab = "Log2(LFQ)",
     main = "Histogram 90 min Log2(LFQ) COVID19 Sample3")

##### Meta Data ######

coagulation_genes = meta$Gene.names

uniprot_genes = vector()
for (i in 1:length(uniprot_list$Gene.names)){
  genelist = strsplit(uniprot_list$Gene.names, " ")
  gene = genelist[[i]][1]
  uniprot_genes = append(uniprot_genes,gene)
}

uniprot_genes[which(coagulation_genes %in% uniprot_genes)]

uniprot_genes[which(coagulation_genes %in% uniprot_genes)]
#### Peptide Data ######
df <- subset_proteins_CoagCascade
df_peptide <- plasma_peptide_coagulationCascade
gene = coagulation_genes[1]



for (i in 1:length(coagulation_genes)) {
  gene = coagulation_genes[i]
  data <- df[grep(gene, df$Fasta.headers),]
  peptides <- df_peptide[grep(gene, df_peptide$Gene.names),]
  
  write.csv(peptides, paste0("P:/All_20200428_COVID_plasma_multiomics/Proteomics/90 min single-shot reps for Anji/CoagulationCascade_peptide_data/",
                             gene,"_peptides.csv"))
}

boxplot(t(clean_df_log2))

pheatmap(clean_df_log2)
