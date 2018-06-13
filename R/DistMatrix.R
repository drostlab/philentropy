DistMatrix <- function(dists, func, testNA, unit){
  
  FUN <- match.fun(func)
  n_cols = ncol(dists)
  dist_value = vector("numeric", 1)
  dist_matrix <- matrix(NA_real_, n_cols,n_cols)
  
  for (i in seq_len(n_cols)) {
    for (j in seq_len(n_cols)) {
      if (is.na(dist_matrix[i,j])) {
        dist_value = FUN(dists[ , i], dists[ , j], testNA, unit)
        dist_matrix[i,j] = dist_value
        dist_matrix[j,i] = dist_value
      }
    }
  } 
  
  return(dist_matrix)
}
