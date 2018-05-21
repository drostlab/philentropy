context("Test implementation of manhattan distance ...")

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



test_that("distance(method = 'manhattan') computes the correct distance value.",
          {
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "manhattan")),
                               sum(abs((P) - (Q))))
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "manhattan")),
                               as.vector(stats::dist(base::rbind(P, Q), method = "manhattan")))
                  
                  expect_equal(as.vector(philentropy::distance(rbind(V, W), method = "manhattan")),
                               sum(abs((V) - (W))))
                  
                  expect_equal(as.vector(philentropy::distance(rbind(V, W), method = "manhattan")),
                               as.vector(stats::dist(base::rbind(V, W), method = "manhattan")))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "manhattan")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               as.vector(dist(distMat, method = "manhattan")))
          })


test_that("Correct manhattan distance is computed when vectors contain 0 values ...", {
        P1 <- c(1,0)
        P2 <- c(0.5, 0.5)
        Q1 <- c(0.5,0.5)
        Q2 <- c(1,0)
        
        expect_equal(as.vector(philentropy::distance(rbind(P1, Q1), method = "manhattan")),
                     as.vector(stats::dist(base::rbind(P1, Q1), method = "manhattan")))
        
        expect_equal(as.vector(philentropy::distance(rbind(P2, Q2), method = "manhattan")),
                     as.vector(stats::dist(base::rbind(P2, Q2), method = "manhattan")))

})


