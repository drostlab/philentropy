context("Test implementation of lorentzian distance ...")

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

test_lorentzian_dist <- function(P, Q) {
  sum(log(1 + abs((P) - (Q))))
}

test_that("distance(method = 'lorentzian') computes the correct distance value using unit = log.",
          {
            expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "lorentzian")),
                         test_lorentzian_dist(P, Q))
            
            # test correct computation of distance matrix
            distMat <-
              rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
            dist.vals <-
              distance(distMat, method = "lorentzian")
            
            expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                         test_dist_matrix(distMat, FUN = test_lorentzian_dist))
            
          })

test_that(
  "distance(method = 'lorentzian') computes the correct distance value using unit = log2.",
  {
    expect_equal(as.vector(
      philentropy::distance(rbind(P, Q), method = "lorentzian",
                            unit = "log2")
    ), sum(log2(1 + abs((
      P
    ) - (
      Q
    )))))
    
  }
)

test_that(
  "distance(method = 'lorentzian') computes the correct distance value using unit = log10.",
  {
    expect_equal(as.vector(
      philentropy::distance(rbind(P, Q), method = "lorentzian",
                            unit = "log10")
    ), sum(log10(1 + abs((
      P
    ) - (
      Q
    )))))
    
    
  }
)


test_that("Correct lorentzian distance is computed when vectors contain 0 values ...",
          {
            P1 <- c(1, 0)
            P2 <- c(0.5, 0.5)
            Q1 <- c(0.5, 0.5)
            Q2 <- c(1, 0)
            
            distMat <-
              rbind(P1, Q1, P2, Q2)
            dist.vals <-
              distance(distMat, method = "lorentzian")
            
            expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                         test_dist_matrix(distMat, FUN = test_lorentzian_dist))
            
          })
