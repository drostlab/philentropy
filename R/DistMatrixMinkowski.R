DistMatrixMinkowski <- function(dists, p, testNA){
  
  n_cols = ncol(dists)
  dist_value = vector("numeric", 1)
  dist_matrix <- matrix(NA_real_, n_cols,n_cols)
  
  for (i in seq_len(n_cols)) {
    for (j in seq_len(n_cols)) {
      if (is.na(dist_matrix[i,j])) {
        dist_value = minkowski(dists[ , i], dists[ , j], p, testNA)
        dist_matrix[i,j] = dist_value
        dist_matrix[j,i] = dist_value
      }
    }
  } 
  
  return(dist_matrix)
}
