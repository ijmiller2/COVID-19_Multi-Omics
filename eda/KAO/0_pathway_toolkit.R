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


