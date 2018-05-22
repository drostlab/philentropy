context("Test implementation of wavehedges distance ...")


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
      res.dist.matrix[i, j] <- dist.fun(x[i, ], x[j, ])
    }
  }
  return(res.dist.matrix[lower.tri(res.dist.matrix, diag = FALSE)])
}



test_that("distance(method = 'wavehedges') computes the correct distance value.",
          {
            test_wh_dist <- function(P, Q) {
              sum(abs(P - Q) / apply(base::rbind(P, Q), 2, max))
            }
            
            expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "wavehedges")),
                         test_wh_dist(P, Q))
            
            # test correct computation of distance matrix
            distMat <-
              rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
            dist.vals <-
              distance(distMat, method = "wavehedges")
            
            expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                         test_dist_matrix(distMat, FUN = test_wh_dist))
          })



test_that(
  "distance(method = 'wavehedges') computes the correct distance value in case the input probability vectors store 0 values at the same position causing 0/0 computations.",
  {
    A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
    B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
    
    wh <- function(x, y) {
      dist <- vector(mode = "numeric", length = 1)
      dist <- 0
      
      for (i in 1:length(x)) {
        if ((abs(x[i] - y[i]) == 0) | ((max(c(
          x[i], y[i]
        ))) == 0)) {
          dist = dist
        } else {
          dist = dist + (abs(x[i] - y[i]) / max(c(x[i], y[i])))
        }
      }
      return(dist)
    }
    
    expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "wavehedges")), wh(A, B))
    
  }
)


