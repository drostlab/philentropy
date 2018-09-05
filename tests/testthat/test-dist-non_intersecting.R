context("Test implementation of non-intersection distance ...")

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


test_that("distance(method = 'non-intersection') computes the correct distance value.",
          {
            expect_equal(as.vector(
              philentropy::distance(rbind(P, Q), method = "non-intersection")
            ),
            1 - sum(apply(base::rbind(P, Q), 2, min)))
            
          })

test_that("Correct non-intersection distance is computed when vectors contain 0 values ...",
          {
            P <- c(1, 0)
            Q <- c(0.5, 0.5)

            expect_equal(as.vector(
              philentropy::distance(rbind(P, Q), method = "non-intersection")
            ),
            1 - sum(apply(base::rbind(P, Q), 2, min)))
            
          })
            

