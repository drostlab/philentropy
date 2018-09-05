context("Test implementation of euclidean distance ...")


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

test_that(
        "'euclidien' (or any other not implemented string) is caught when wrong input string for method is entered",
        {
                distMat <- rbind(rep(0.2, 5), rep(0.1, 10))
                expect_error(
                        distance(distMat, method = "euclidien"),
                        "Method 'euclidien' is not implemented in this function. Please consult getDistMethods()."
                )
        }
)

test_that("Only numeric values are passed to distance()", {
        distMat <- rbind(rep("A", 10), rep("B", 10))
        expect_error(
                distance(distMat, method = "euclidean"),
                paste0(
                        "Your input ",
                        class(distMat),
                        " stores non-numeric values. Non numeric values cannot be used to compute distances.."
                )
        )
})


test_that("Only choose from units: log, log2, or log10", {
        distMat <- rbind(rep(0.2, 5), rep(0.1, 10))
        expect_error(
                distance(distMat, method = "euclidean", unit = "log5"),
                "You can only choose units: log, log2, or log10."
        )
})

test_that("distance(method = 'euclidean') computes the correct distance value.",
          {
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "euclidean")),
                               sqrt(sum(abs((
                                       P
                               ) - (
                                       Q
                               )) ^ 2)))
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "euclidean")),
                               as.vector(stats::dist(base::rbind(P, Q), method = "euclidean")))
                  
                  expect_equal(as.vector(philentropy::distance(rbind(V, W), method = "euclidean")),
                               as.vector(stats::dist(base::rbind(V, W), method = "euclidean")))
                  
                  expect_equal(as.vector(philentropy::distance(rbind(V, W), method = "euclidean")),
                               sqrt(sum(abs((
                                       V
                               ) - (
                                       W
                               )) ^ 2)))
                  
                  
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "euclidean")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               as.vector(dist(distMat)))
                  
                  #expect_error(philentropy::distance(1:10, 20:29, method = "euclidean"))
          })


test_that("Correct euclidean distance is computed when vectors contain 0 values ...", {
        P1 <- c(1,0)
        P2 <- c(0.5, 0.5)
        Q1 <- c(0.5,0.5)
        Q2 <- c(1,0)
        
        expect_equal(as.vector(philentropy::distance(rbind(P1, Q1), method = "euclidean")),
                     as.vector(stats::dist(base::rbind(P1, Q1), method = "euclidean")))
        
        expect_equal(as.vector(philentropy::distance(rbind(P2, Q2), method = "euclidean")),
                     as.vector(stats::dist(base::rbind(P2, Q2), method = "euclidean")))
        
        distMat <- rbind(P1, Q1, P2, Q2)
        dist.vals <-
                distance(distMat, method = "euclidean")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     as.vector(dist(distMat)))
})






