context("Test implementation of minkowski distance ...")


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


test_that("distance(method = 'minkowski') computes the correct distance value.",
          {
                  expect_equal(as.vector(
                          philentropy::distance(rbind(P, Q), method = "minkowski", p = 4)
                  ),
                  (sum(abs((
                          P
                  ) - (
                          Q
                  )) ^ 4)) ^ 0.25)
                  
                  expect_equal(as.vector(
                          philentropy::distance(rbind(P, Q), method = "minkowski", p = 4)
                  ),
                  as.vector(stats::dist(
                          base::rbind(P, Q), method = "minkowski", p = 4
                  )))
                  
                  expect_error(
                          as.vector(philentropy::distance(rbind(P, Q), method = "minkowski")),
                          "Please specify p for the Minkowski distance."
                  )
                  
                  expect_equal(as.vector(
                          philentropy::distance(rbind(V, W), method = "minkowski", p = 4)
                  ), (sum(abs((
                          V
                  ) - (
                          W
                  )) ^ 4)) ^ 0.25)
                  
                  expect_equal(as.vector(
                          philentropy::distance(rbind(V, W), method = "minkowski", p = 4)
                  ),
                  as.vector(stats::dist(
                          base::rbind(V, W), method = "minkowski", p = 4
                  )))
                  
                  expect_error(
                          as.vector(philentropy::distance(rbind(V, W), method = "minkowski")),
                          "Please specify p for the Minkowski distance."
                  )
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "minkowski", p = 4)
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               as.vector(dist(
                                       distMat, method = "minkowski", p = 4
                               )))
                  
          })



test_that("Correct minkowski distance is computed when vectors contain 0 values ...", {
        P1 <- c(1,0)
        P2 <- c(0.5, 0.5)
        Q1 <- c(0.5,0.5)
        Q2 <- c(1,0)
        
        expect_equal(as.vector(philentropy::distance(rbind(P1, Q1), method = "minkowski", p = 4)),
                     as.vector(stats::dist(base::rbind(P1, Q1), method = "minkowski", p = 4)))
        
        expect_equal(as.vector(philentropy::distance(rbind(P2, Q2), method = "minkowski", p = 4)),
                     as.vector(stats::dist(base::rbind(P2, Q2), method = "minkowski", p = 4)))
        
})
