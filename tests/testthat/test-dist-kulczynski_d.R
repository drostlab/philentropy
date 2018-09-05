context("Test implementation of kulczynski_d distance ...")


P <- 1:10 / sum(1:10)
Q <- 20:29 / sum(20:29)
V <- -10:10
W <- -20:0

# function to test distance matrix functionality
# for different distance measures
test_dist_matrix <- function(x, FUN) {
  dist.fun <- match.fun(FUN)
  res.dist.matrix <- matrix(NA_real_, nrow(x), nrow(x))
  
  for (i in 1:nrow(x)) {
    for (j in 1:nrow(x)) {
      res.dist.matrix[i, j] <- dist.fun(x[i,], x[j,])
    }
  }
  return(res.dist.matrix[lower.tri(res.dist.matrix, diag = FALSE)])
}


test_kd_dist <- function(P, Q) {
  #sum(abs((P) - (Q))) / sum(apply(rbind(P, Q), 2, min))
   diff       = 0.0;
   dist1      = 0.0;
   dist2      = 0.0;
   min_point  = 0.0;
  
  for (i in seq_len(length(P))) {
    diff      = abs(P[i] - Q[i]);
    if (P[i] <= Q[i]){
      min_point = P[i];
    } else {
      min_point = Q[i];
    }
    dist1 = diff + dist1;
    if (min_point == 0.0){
      dist2 = dist2 + 0.00001;
    } else {
      dist2 = dist2 + min_point;
    }     
  }
   
   if (dist2 == 0.0) {
     return (NaN);
   } else {
     return (dist1/dist2);         
   }    
}


test_that("distance(method = 'kulczynski_d') computes the correct distance value.",
          {
            expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "kulczynski_d")),
                         test_kd_dist(P, Q))
            
            # test correct computation of distance matrix
            distMat <-
              rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
            dist.vals <-
              distance(distMat, method = "kulczynski_d")
            
            expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                         test_dist_matrix(distMat, FUN = test_kd_dist))
            
          })


test_that("Correct kulczynski_d distance is computed when vectors contain 0 values ...",
          {
            P1 <- c(1, 0)
            P2 <- c(0.5, 0.5)
            Q1 <- c(0.5, 0.5)
            Q2 <- c(1, 0)
            
            distMat <-
              rbind(P1, Q1, P2, Q2)
            dist.vals <-
              distance(distMat, method = "kulczynski_d")
            
            expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                         test_dist_matrix(distMat, FUN = test_kd_dist))
            
          })
