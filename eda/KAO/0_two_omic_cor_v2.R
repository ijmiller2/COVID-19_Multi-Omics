######### Function for correlations #########
two_omic_cor <- function(x, y, method = "pearson") {
  # x is a matrix that will be presented in rows
  # y is a matrix that will be presented in columns
  # x and y must have equal column length 
  if(ncol(x) != ncol(x)) stop("Number of columns are not equal")
  
  x_y_matrix <- matrix(NA, nrow = nrow(x), ncol = nrow(y), dimnames = list(row.names(x), row.names(y)))
  x_y_pvalues <- matrix(NA, nrow = nrow(x), ncol = nrow(y), dimnames = list(row.names(x), row.names(y)))
  
  for(i in 1:nrow(x)){
    for (k in 1:nrow(y)){
      cortest <- cor.test(x[i,], y[k,], method = method)
      x_y_matrix[i,k] <- cortest$estimate[[1]]
      x_y_pvalues[i,k] <- cortest$p.value
    }
  }  
  list(cor = x_y_matrix, pvalue = x_y_pvalues, adjusted_pvalue = matrix(p.adjust(x_y_pvalues), nrow = nrow(x), ncol = nrow(y)))
}


all_omic_cor <- function(x) {
  cor_matrix_long <- NULL

  for(i in 1:nrow(x)){
    for (k in 1:nrow(x)){
     if(i > k){
       cortest <- cor.test(x[i,], x[k,])
       cor_matrix_long <- rbind(cor_matrix_long, c(i, k, cortest$estimate[[1]], cortest$p.value))
     }
      
    }
  }
  cor_matrix_long <-as.data.frame(cor_matrix_long)
  cor_matrix_long[,1] <- row.names(x)[cor_matrix_long[,1]]
  
  cor_matrix_long[,2] <- row.names(x)[cor_matrix_long[,2]]
  names(cor_matrix_long) <- c("compound1", "compound2", "pearson", "pvalue")
  cor_matrix_long
  
}


two_omic_cor_fast <- function(x, y, method = "pearson") {
  # x is a matrix that will be presented in rows
  # y is a matrix that will be presented in columns
  # x and y must have equal column length 
  if(ncol(x) != ncol(x)) stop("Number of columns are not equal")
  apply(x, 1, function(x_i) apply(y, 1, function(y_i) cor(x_i, y_i, method = method)))

}


