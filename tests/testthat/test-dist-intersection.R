context("Test implementation of intersection distance ...")


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

test_intersection_dist <- function(P, Q) {
  sum(apply(base::rbind(P, Q), 2, min))
}

test_that("distance(method = 'intersection') computes the correct distance value.",
          {
            
            
            expect_equal(as.vector(
              philentropy::distance(rbind(P, Q), method = "intersection")
            ),
            test_intersection_dist(P, Q))
            
            # test correct computation of distance matrix
            distMat <-
              rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
            dist.vals <-
              distance(distMat, method = "intersection")
            
            expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                         test_dist_matrix(distMat, FUN = test_intersection_dist))
            
          })



test_that("Correct intersection distance is computed when vectors contain 0 values ...",
          {
            P1 <- c(1, 0)
            P2 <- c(0.5, 0.5)
            Q1 <- c(0.5, 0.5)
            Q2 <- c(1, 0)
            
            distMat <-
              rbind(P1, Q1, P2, Q2)
            dist.vals <-
              distance(distMat, method = "intersection")
            
            expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                         test_dist_matrix(distMat, FUN = test_intersection_dist))
            
          })

