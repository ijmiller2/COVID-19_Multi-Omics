## Use Katie's code for GO Analysis on proteins measured in proteomics
## samples of COVID19

library(DBI)
library(RSQLite)
library(pheatmap)
library(plyr)
library(RColorBrewer)
library(DelayedMatrixStats)
library(randomcoloR)
library(cluster)

#### Upload Data ####
GO_df <- read.csv("H:/Projects/COVID19/Proteomics/Files/uniprot_GO_BioProcess_CellularComp_MolecularFun.csv", 
                  header = TRUE, sep = ",", stringsAsFactors = FALSE)

GO_df_biologicalProcess <- GO_df[,c(1,4)] 
colnames(GO_df_biologicalProcess) <- c("UniprotProtein_First","GO_biologicalProcess")


con <- dbConnect(SQLite(), dbname = "P:/All_20200428_COVID_plasma_multiomics/SQLite Database/Covid-19 Study DB.sqlite")
proteins_meta <- dbGetQuery(con, "SELECT standardized_name, metadata.biomolecule_id, metadata_type, metadata_value
                               FROM metadata
                               INNER JOIN biomolecules ON biomolecules.biomolecule_id = metadata.biomolecule_id
                               WHERE biomolecules.omics_id = 1
                               AND biomolecules.keep = 1
                               ")
dbDisconnect(con)

##### Protein meta formattting #####

# Isolate only the first protein in protein group for Majority.protein.IDs and Protein.IDs
majorityProtein_isoform <- lapply(strsplit(proteins_meta$standardized_name, ";"), '[', 1)
majorityProtein <- sub("-.*", "", majorityProtein_isoform)
proteins_meta$majority.protein <- vapply(majorityProtein, paste, collapse = ", ", character(1L))
write.csv(proteins_meta, "H:/Projects/COVID19/Proteomics/Files/Proteins_meta_MajorityProtein.csv")

# subset metadata with "gene_name"
Go_proteins_meta <- proteins_meta[which(proteins_meta$metadata_type == "gene_name"),]

#which(Go_proteins_meta$majority.protein %in% GO_df_biologicalProcess$UniprotProtein_First)

# For every uniprot entry in GO df check to see if there was a GO term associated. 
protein <- vector()
missing_protein = integer()
for(i in 1:nrow(GO_df_biologicalProcess)){
  if(GO_df_biologicalProcess$GO_biologicalProcess[i] != ""){
    protein <- append(protein, GO_df_biologicalProcess$UniprotProtein_First[i])
  }else{print(paste0("i: ", i, " ", GO_df_biologicalProcess$UniprotProtein_First[i]))
        missing_protein = append(missing_protein,GO_df_biologicalProcess$UniprotProtein_First[i])}
}

which(!protein %in% Go_proteins_meta$majority.protein)

biomol_id <- vector()
for(i in 1:length(protein)){
  extract_protein <- protein[i]
  print(paste0("i: ", i, " ", extract_protein))
  #id <- Go_proteins_meta[grep(extract_protein, Go_proteins_meta$majority.protein),2]
  id <- Go_proteins_meta[which(Go_proteins_meta$majority.protein %in% extract_protein),2]
  biomol_id <- append(biomol_id, id)
}


Go_proteins_notmissing <- Go_proteins_meta[which(Go_proteins_meta$biomolecule_id %in% biomol_id),]
GO_df_biologicalProcess_complete <- GO_df_biologicalProcess[which(protein %in% GO_df_biologicalProcess$UniprotProtein_First),]

reference_index <- GO_df_biologicalProcess_complete$UniprotProtein_First

id_index <- GO_df_biologicalProcess_complete$GO_biologicalProcess

############ Creating pathway analysis tool kit #############

############ Function for making reference sets #############


make_reference_sets <- function(reference_index, id_index, split_str = "; "){
  all_reference <- NULL
  for(i in 1:length(reference_index)){
    all_reference <- append(all_reference, strsplit(reference_index[i], split_str)[[1]])
  }
  
  unique_reference <- as.list(unique(all_reference), stringsAsFactors = F)
  
  reference_sets <- lapply(unique_reference, function(x) id_index[grep(x[1], reference_index, fixed= T)])
  names(reference_sets) <- unique_reference
  
  reference_sets
}

reference_sets <- make_reference_sets(reference_index, id_index, split_str = ";")

############ Function for testing enrichment ############## 

enrichment <- function(set, reference_sets, background){
  # output p_value
  # output fdr
  nset <- length(set)
  set <- as.character(set)
  nbackground <- length(background)
  background <- as.character(background)  
  output <- data.frame(reference = names(reference_sets), 
                       pvalue = rep(1, length(names(reference_sets))), 
                       fdr_pvalue = rep(NA, length(names(reference_sets))),
                       stringsAsFactors = F)
  for (i in 1:nrow(output)){
    hits <- length(intersect(set, reference_sets[[output[i,1]]]))
    hitsBackground <- length(intersect(background, reference_sets[[output[i,1]]]))
    if (hits > 0){
      output$pvalue[i] <- phyper(hits-1, hitsBackground, length(background) - hitsBackground, length(set), lower.tail = F)
    } 
    if (length(reference_sets[[output[i,1]]]) == 1){
      output$pvalue[i] <- NA
    }
  }
  output$fdr_pvalue <- p.adjust(output$pvalue, method = "BH")
  output
} 

set <- Deane_proteins$UniProt.ID
background <- reference_index
COVID_enrichment <- enrichment(set, reference_sets, background)

