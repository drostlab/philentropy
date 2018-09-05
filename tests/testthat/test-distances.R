#  Part of the philentropy package
#
#  Copyright (C) 2015 Hajk-Georg Drost
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

context("Test implementation of distance measures...")



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

context("Test implementation of motyka distance ...")

test_that("distance(method = 'motyka') computes the correct distance value.",
          {
                  test_motyka_dist <- function(P, Q) {
                          sum(apply(base::rbind(P, Q), 2, max)) / sum(P + Q)
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "motyka")),
                               test_motyka_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <- distance(distMat, method = "motyka")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_motyka_dist))
          })


context("Test implementation of kulczynski_s distance ...")

test_that("distance(method = 'kulczynski_s') computes the correct distance value.",
          {
                  test_ks_dist <- function(P, Q) {
                          1 / (sum(abs((P) - (Q))) / sum(apply(rbind(
                                  P, Q
                          ), 2, min)))
                  }
                  
                  expect_equal(as.vector(
                          philentropy::distance(rbind(P, Q), method = "kulczynski_s")
                  ),
                  test_ks_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "kulczynski_s")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_ks_dist))
                  
          })

context("Test implementation of tanimoto distance ...")

test_that("distance(method = 'tanimoto') computes the correct distance value.",
          {
                  test_tanimoto_dist <- function(P, Q) {
                          sum(abs((P) - (Q))) / sum(apply(rbind(P, Q), 2, max))
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "tanimoto")),
                               test_tanimoto_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "tanimoto")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_tanimoto_dist))
          })


context("Test implementation of ruzicka distance ...")

test_that("distance(method = 'ruzicka') computes the correct distance value.",
          {
                  test_ruzicka_dist <- function(P, Q) {
                          1 - (sum(abs((P) - (Q))) / sum(apply(rbind(
                                  P, Q
                          ), 2, max)))
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "ruzicka")),
                               test_ruzicka_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <- distance(distMat, method = "ruzicka")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_ruzicka_dist))
                  
          })


context("Test implementation of inner_product distance ...")

test_that("distance(method = 'inner_product') computes the correct distance value.",
          {
                  test_innerproduct_dist <- function(P, Q) {
                          sum (P * Q)
                  }
                  
                  expect_equal(as.vector(
                          philentropy::distance(rbind(P, Q), method = "inner_product")
                  ),
                  test_innerproduct_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "inner_product")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_innerproduct_dist))
          })


context("Test implementation of harmonic_mean distance ...")

test_that("distance(method = 'harmonic_mean') computes the correct distance value.",
          {
                  test_harmonic_mean_dist <- function(P, Q) {
                          2 * sum ((P * Q) / (P + Q))
                  }
                  
                  expect_equal(as.vector(
                          philentropy::distance(rbind(P, Q), method = "harmonic_mean")
                  ),
                  test_harmonic_mean_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "harmonic_mean")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_harmonic_mean_dist))
                  
          })


test_that(
        "distance(method = 'harmonic_mean') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0/0 computation
        
        .",
        {
                A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                
                hm <- function(x, y) {
                        dist <- vector(mode = "numeric", length = 1)
                        dist <- 0
                        
                        for (i in 1:length(x)) {
                                if (((x[i] * y[i]) == 0) & ((x[i] + y[i]) == 0)) {
                                        dist = dist
                                } else {
                                        dist = dist + ((x[i] * y[i]) / (x[i] + y[i]))
                                }
                                
                        }
                        
                        return(2 * dist)
                }
                
                expect_equal(as.vector(
                        philentropy::distance(rbind(A, B), method = "harmonic_mean")
                ), hm(A, B))
        }
)


context("Test implementation of cosine distance ...")

test_that("distance(method = 'cosine') computes the correct distance value.",
          {
                  test_cosine_dist <- function(P, Q) {
                          ((sum ((P) * (Q))) / (sqrt(sum((P) ^ 2
                          )) * sqrt(sum((Q) ^ 2
                          ))))
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "cosine")),
                               test_cosine_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <- distance(distMat, method = "cosine")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_cosine_dist))
                  
          })


context("Test implementation of hassebrook distance ...")

test_that("distance(method = 'hassebrook') computes the correct distance value.",
          {
                  test_hassebrook_dist <- function(P, Q) {
                          ((sum ((P) * (Q))) / (sum((P) ^ 2) + sum((Q) ^ 2) - ((
                                  sum ((P) * (Q))
                          ))))
                  }
                  
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "hassebrook")),
                               test_hassebrook_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "hassebrook")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_hassebrook_dist))
                  
          })


context("Test implementation of jaccard distance ...")

test_that("distance(method = 'jaccard') computes the correct distance value.",
          {
                  test_jaccard_dist <- function(P, Q) {
                          1 - ((sum ((P) * (Q))) / (sum((P) ^ 2) + sum((Q) ^ 2) - ((
                                  sum ((P) * (Q))
                          ))))
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "jaccard")),
                               test_jaccard_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <- distance(distMat, method = "jaccard")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_jaccard_dist))
          })


context("Test implementation of dice distance ...")

test_that("distance(method = 'dice') computes the correct distance value.",
          {
                  test_dice_dist <- function(P, Q) {
                          1 - (2 * (sum ((P) * (Q))) / (sum((P) ^ 2) + sum((Q) ^ 2)))
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "dice")),
                               test_dice_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <- distance(distMat, method = "dice")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_dice_dist))
                  
          })


context("Test implementation of fidelity distance ...")

test_that("distance(method = 'fidelity') computes the correct distance value.",
          {
                  test_fidelity_dist <- function(P, Q) {
                          sum(sqrt(P * Q))
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "fidelity")),
                               test_fidelity_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "fidelity")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_fidelity_dist))
                  
          })


context("Test implementation of bhattacharyya distance ...")


test_that(
        "distance(method = 'bhattacharyya') computes the correct distance value using unit = log.",
        {
                test_b_dist <- function(P, Q) {
                        -log(sum(sqrt(P * Q)))
                }
                
                expect_equal(as.vector(
                        philentropy::distance(rbind(P, Q), method = "bhattacharyya")
                ),
                test_b_dist(P, Q))
                
                # test correct computation of distance matrix
                distMat <-
                        rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                dist.vals <-
                        distance(distMat, method = "bhattacharyya")
                
                expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                             test_dist_matrix(distMat, FUN = test_b_dist))
                
        }
)



test_that(
        "distance(method = 'bhattacharyya') computes the correct distance value using unit = log2.",
        {
                expect_equal(as.vector(
                        philentropy::distance(rbind(P, Q), method = "bhattacharyya", unit = "log2")
                ), -log2(sum(sqrt(P * Q))))
                
        }
)

test_that(
        "distance(method = 'bhattacharyya') computes the correct distance value using unit = log10.",
        {
                expect_equal(as.vector(
                        philentropy::distance(rbind(P, Q), method = "bhattacharyya", unit = "log10")
                ), -log10(sum(sqrt(P * Q))))
                
        }
)


context("Test implementation of hellinger distance ...")

test_that("distance(method = 'hellinger') computes the correct distance value.",
          {
                  test_hellinger_dist <- function(P, Q) {
                          2L * sqrt(1L - sum(sqrt(P * Q)))
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "hellinger")),
                               test_hellinger_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5) / sum(c(5, 1, 7, 9, 5)))
                  dist.vals <-
                          distance(distMat, method = "hellinger")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_hellinger_dist))
                  
          })



context("Test implementation of matusita distance ...")

test_that("distance(method = 'matusita') computes the correct distance value.",
          {
                  test_matusita_dist <- function(P, Q) {
                          sqrt(sum((sqrt(P) - sqrt(Q)) ^ 2))
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "matusita")),
                               test_matusita_dist(P, Q))
                  
                  
                  # test correct computation of distance matrix
                  A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                  B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                  C <- c(0, 0.25, 0, 0.25, 0.25, 0.25)
                  distMat <- rbind(A, B, C)
                  dist.vals <-
                          distance(distMat, method = "matusita")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_matusita_dist))
          })



context("Test implementation of squared_chord distance ...")

test_that("distance(method = 'squared_chord') computes the correct distance value.",
          {
                  test_sqchord_dist <- function(P, Q) {
                          sum((sqrt(P) - sqrt(Q)) ^ 2)
                  }
                  
                  expect_equal(as.vector(
                          philentropy::distance(rbind(P, Q), method = "squared_chord")
                  ),
                  test_sqchord_dist(P, Q))
                  
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "squared_chord")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_sqchord_dist))
          })



context("Test implementation of squared_euclidean distance ...")

test_that("distance(method = 'squared_euclidean') computes the correct distance value.",
          {
                  test_sqeuclidean_dist <- function(P, Q) {
                          sum(((P) - (Q)) ^ 2)
                  }
                  
                  expect_equal(as.vector(
                          philentropy::distance(rbind(P, Q), method = "squared_euclidean")
                  ),
                  test_sqeuclidean_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "squared_euclidean")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_sqeuclidean_dist))
                  
          })



context("Test implementation of pearson distance ...")


test_that("distance(method = 'pearson') computes the correct distance value.",
          {
                  test_pearson_dist <- function(P, Q) {
                          sum((P - Q) ^ 2 / Q)
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "pearson")),
                               test_pearson_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                  B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                  C <- c(0, 0.25, 0, 0.25, 0.25, 0.25)
                  distMat <- rbind(A, B, C)
                  dist.vals <- distance(distMat, method = "pearson")
                  
                  #         expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                  #                      test_dist_matrix(distMat, FUN = test_pearson_dist))
                  
          })


context("Test implementation of neyman distance ...")


test_that("distance(method = 'neyman') computes the correct distance value.",
          {
                  test_neyman_dist <- function(P, Q) {
                          sum(((P - Q) ^ 2) / P)
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "neyman")),
                               test_neyman_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                  B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                  C <- c(0, 0.25, 0, 0.25, 0.25, 0.25)
                  distMat <- rbind(A, B, C)
                  dist.vals <- distance(distMat, method = "neyman")
                  
                  #         expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                  #                      test_dist_matrix(distMat, FUN = test_neyman_dist))
          })


context("Test implementation of squared_chi distance ...")


test_that("distance(method = 'squared_chi') computes the correct distance value.",
          {
                  test_sqchi_dist <- function(P, Q) {
                          sum(((P) - (Q)) ^ 2 / ((P) + (Q)))
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "squared_chi")),
                               test_sqchi_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "squared_chi")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_sqchi_dist))
          })



test_that(
        "distance(method = 'squared_chi') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0/0 computation
        .",
        {
                A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                
                sqchisq <- function(x, y) {
                        dist <- vector(mode = "numeric", length = 1)
                        dist <- 0
                        
                        for (i in 1:length(x)) {
                                if (((x[i] - y[i]) ^ 2 == 0) & ((x[i] + y[i]) == 0)) {
                                        dist = dist
                                } else {
                                        dist = dist + ((x[i] - y[i]) ^ 2 / (x[i] + y[i]))
                                }
                                
                        }
                        
                        return(dist)
                }
                
                expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "squared_chi")), sqchisq(A, B))
                
        }
)


context("Test implementation of prob_symm distance ...")

test_that("distance(method = 'prob_symm') computes the correct distance value.",
          {
                  test_probsymm_dist <- function(P, Q) {
                          2 * sum(((P) - (Q)) ^ 2 / ((P) + (Q)))
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "prob_symm")),
                               test_probsymm_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "prob_symm")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_probsymm_dist))
          })




test_that(
        "distance(method = 'prob_symm') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0/0 computation
        .",
        {
                A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                
                probsymmchisq <- function(x, y) {
                        dist <- vector(mode = "numeric", length = 1)
                        dist <- 0
                        
                        for (i in 1:length(x)) {
                                if (((x[i] - y[i]) ^ 2 == 0) & ((x[i] + y[i]) == 0)) {
                                        dist = dist
                                } else {
                                        dist = dist + ((x[i] - y[i]) ^ 2 / (x[i] + y[i]))
                                }
                                
                        }
                        
                        return(2 * dist)
                }
                
                expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "prob_symm")),
                             probsymmchisq(A, B))
                
        }
)


context("Test implementation of divergence distance ...")

test_that("distance(method = 'divergence') computes the correct distance value.",
          {
                  test_divergence_dist <- function(P, Q) {
                          2 * sum(((P) - (Q)) ^ 2 / ((P) + (Q)) ^ 2)
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "divergence")),
                               test_divergence_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "divergence")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_divergence_dist))
          })


test_that(
        "distance(method = 'divergence') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0/0 computation
        .",
        {
                A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                
                div <- function(x, y) {
                        dist <- vector(mode = "numeric", length = 1)
                        dist <- 0
                        
                        for (i in 1:length(x)) {
                                if (((x[i] - y[i]) ^ 2 == 0) & ((x[i] + y[i]) ^ 2 == 0)) {
                                        dist = dist
                                } else {
                                        dist = dist + ((x[i] - y[i]) ^ 2 / (x[i] + y[i]) ^ 2)
                                }
                                
                        }
                        
                        return(2 * dist)
                }
                
                expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "divergence")), div(A, B))
        }
)


context("Test implementation of clark distance ...")


test_that("distance(method = 'clark') computes the correct distance value.",
          {
                  test_clark_dist <- function(P, Q) {
                          sqrt(sum((abs((P) - (Q)
                          ) / ((P) + (Q)
                          )) ^ 2))
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "clark")),
                               test_clark_dist(P, Q))
                  
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <- distance(distMat, method = "clark")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_clark_dist))
          })


test_that(
        "distance(method = 'clark') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0/0 computation
        .",
        {
                A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                
                clark <- function(x, y) {
                        dist <- vector(mode = "numeric", length = 1)
                        dist <- 0
                        
                        for (i in 1:length(x)) {
                                if ((abs(x[i] - y[i]) == 0) & ((x[i] + y[i]) == 0)) {
                                        dist = dist
                                } else {
                                        dist = dist + (abs(x[i] - y[i]) / (x[i] + y[i])) ^ 2
                                }
                                
                        }
                        
                        return(sqrt(dist))
                }
                
                expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "clark")), clark(A, B))
                
        }
)



context("Test implementation of additive_symm distance ...")

test_that("distance(method = 'additive_symm') computes the correct distance value.",
          {
                  test_addsymm_dist <- function(P, Q) {
                          sum((((P) - (Q)) ^ 2 * ((P) + (Q))) / ((P) * (Q)))
                  }
                  
                  expect_equal(as.vector(
                          philentropy::distance(rbind(P, Q), method = "additive_symm")
                  ),
                  test_addsymm_dist(P, Q))
                  
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "additive_symm")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_addsymm_dist))
          })


test_that(
        "distance(method = 'additive_symm') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0/0 computation
        .",
        {
                A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                
                add <- function(x, y) {
                        dist <- vector(mode = "numeric", length = 1)
                        dist <- 0
                        
                        for (i in 1:length(x)) {
                                if (((x[i] + y[i]) == 0) | ((x[i] * y[i]) == 0)) {
                                        dist = dist
                                } else {
                                        dist = dist + ((x[i] - y[i]) ^ 2 * ((x[i] + y[i]) / (x[i] * y[i])))
                                }
                                
                        }
                        
                        return(dist)
                }
                
                expect_equal(as.vector(
                        philentropy::distance(rbind(A, B), method = "additive_symm")
                ), add(A, B))
        }
)


context("Test implementation of kullback-leibler distance ...")


test_that(
        "distance(method = 'kullback-leibler') computes the correct distance value using unit = log.",
        {
                test_KL_dist <- function(P, Q) {
                        dist <- 0
                        for (i in seq_len(length(P))) {
                                
                                if (P[i] == 0 & Q[i] == 0){
                                        dist = dist
                                } else {
                                        if (P[i] == 0) {
                                                dist <- dist + 0
                                        } else {
                                                if (Q[i] == 0) {
                                                        dist <- dist + (P[i] * log(P[i] / 0.00001))
                                                } else {
                                                        dist <- dist + (P[i] * log(P[i] / Q[i]))
                                                }
                                        }
                                }
                                }
                                
                        return(dist)        
                }
                
                expect_equal(as.vector(
                        philentropy::distance(rbind(P, Q), method = "kullback-leibler")
                ),
                test_KL_dist(P, Q))
                
                # test correct computation of distance matrix
                A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                C <- c(0, 0.25, 0, 0.25, 0.25, 0.25)
                distMat <- rbind(A, B, C)
                dist.vals <-
                        distance(distMat, method = "kullback-leibler")
                
                expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                             test_dist_matrix(distMat, FUN = test_KL_dist))
        }
)


test_that(
        "distance(method = 'kullback-leibler') computes the correct distance value using unit = log2.",
        {
                expect_equal(as.vector(
                        philentropy::distance(rbind(P, Q), method = "kullback-leibler", unit = "log2")
                ), sum((P) * log2((P) / (Q))))
                
        }
)

test_that(
        "distance(method = 'kullback-leibler') computes the correct distance valueusing unit = log10.",
        {
                expect_equal(as.vector(
                        philentropy::distance(rbind(P, Q), method = "kullback-leibler", unit = "log10")
                ), sum((P) * log10((P) / (Q))))
                
        }
)

test_that(
        "distance(method = 'kullback-leibler') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0 * log(0) computation
        .",
        {
                A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                
                kl <- function(P, Q) {
                        dist <- 0
                        for (i in seq_len(length(P))) {
                                
                                if (P[i] == 0 & Q[i] == 0){
                                        dist = dist
                                } else {
                                        if (P[i] == 0) {
                                                dist <- dist + 0
                                        } else {
                                                if (Q[i] == 0) {
                                                        dist <- dist + (P[i] * log(P[i] / 0.00001))
                                                } else {
                                                        dist <- dist + (P[i] * log(P[i] / Q[i]))
                                                }
                                        }
                                }
                        }
                        
                        return(dist)        
                }
                
                expect_equal(as.vector(
                        philentropy::distance(rbind(A, B), method = "kullback-leibler")
                ), kl(A, B))
                
        }
)


context("Test implementation of jeffreys distance ...")


test_that("distance(method = 'jeffreys') computes the correct distance value using unit = log.",
          {
                  test_jeffreys_dist <- function(P, Q) {
                          sum(((P) - (Q)) * log((P) / (Q)))
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "jeffreys")),
                               test_jeffreys_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "jeffreys")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_jeffreys_dist))
                  
          })


test_that("distance(method = 'jeffreys') computes the correct distance value using unit = log2.",
          {
                  expect_equal(as.vector(
                          philentropy::distance(rbind(P, Q), method = "jeffreys", unit = "log2")
                  ), sum(((P) - (Q)) * log2((P) / (Q))))
                  
          })

test_that("distance(method = 'jeffreys') computes the correct distance value using unit = log10.",
          {
                  expect_equal(as.vector(
                          philentropy::distance(rbind(P, Q), method = "jeffreys", unit = "log10")
                  ), sum(((P) - (Q)) * log10((P) / (Q))))
                  
          })

test_that(
        "distance(method = 'jeffreys') computes the correct distance value in case 0 values are stored in the input probability vectors -> leading to log(0) computations.",
        {
                A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- rep(1, 6) / 6
                
                jeff <- function(x, y) {
                        if (any((x / y) == 0)) {
                                xy.ratio <- x / y
                                xy.ratio[xy.ratio == 0] <- 0.00001
                                sum((x - y) * log(xy.ratio))
                        } else {
                                sum((x - y) * log(x / y))
                        }
                        
                }
                expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "jeffreys")), jeff(A, B))
        }
)


context("Test implementation of k_divergence distance ...")

test_that(
        "distance(method = 'k_divergence') computes the correct distance value using unit = log.",
        {
                test_kdivergence_dist <- function(P, Q) {
                        sum((P) * log(2 * (P) / ((P) + (Q))))
                }
                
                expect_equal(as.vector(
                        philentropy::distance(rbind(P, Q), method = "k_divergence")
                ),
                test_kdivergence_dist(P, Q))
                
                # test correct computation of distance matrix
                A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                C <- c(0, 0.25, 0, 0.25, 0.25, 0.25)
                
                distMat <- rbind(A, B, C)
                dist.vals <-
                        distance(distMat, method = "k_divergence")
                
                expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                             test_dist_matrix(distMat, FUN = test_kdivergence_dist))
        }
)

test_that(
        "distance(method = 'k_divergence') computes the correct distance value using unit = log2.",
        {
                expect_equal(as.vector(
                        philentropy::distance(rbind(P, Q), method = "k_divergence", unit = "log2")
                ), sum((P) * log2(2 * (P) / ((
                        P
                ) + (
                        Q
                )))))
                
        }
)

test_that(
        "distance(method = 'k_divergence') computes the correct distance value using unit = log10.",
        {
                expect_equal(as.vector(
                        philentropy::distance(rbind(P, Q), method = "k_divergence", unit = "log10")
                ), sum((P) * log10(2 * (P) / ((
                        P
                ) + (
                        Q
                )))))
                
        }
)

test_that(
        "distance(method = 'k_divergence') computes the correct distance value in case 0 values are stored in the input probability vectors -> leading to 0 * log(0) computations.",
        {
                A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                
                kdiv <- function(x, y) {
                        dist <- vector(mode = "numeric", length = 1)
                        dist <- 0
                        
                        for (i in 1:length(x)) {
                                if ((x[i] == 0) & ((y[i]) == 0)) {
                                        dist = dist
                                } else {
                                        dist = dist + (x[i] * log((2 * x[i]) / (x[i] + y[i])))
                                }
                                
                        }
                        
                        return(dist)
                        
                }
                expect_equal(as.vector(
                        philentropy::distance(rbind(A, B), method = "k_divergence")
                ), kdiv(A, B))
        }
)


context("Test implementation of topsoe distance ...")

test_that("distance(method = 'topsoe') computes the correct distance value using unit = log.",
          {
                  test_topsoe_dist <- function(P, Q) {
                          sum(((P) * log(2 * (P) / ((P) + (Q)
                          ))) + ((Q) * log(2 * (Q) / ((P) + (Q)
                          ))))
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "topsoe")),
                               test_topsoe_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <- distance(distMat, method = "topsoe")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_topsoe_dist))
                  
          })


test_that("distance(method = 'topsoe') computes the correct distance value using unit = log2.",
          {
                  expect_equal(as.vector(
                          philentropy::distance(rbind(P, Q), method = "topsoe", unit = "log2")
                  ),
                  sum(((P) * log2(
                          2 * (P) / ((P) + (Q))
                  )) + ((Q) * log2(
                          2 * (Q) / ((P) + (Q))
                  ))))
                  
          })


test_that("distance(method = 'topsoe') computes the correct distance value using unit = log10.",
          {
                  expect_equal(as.vector(
                          philentropy::distance(rbind(P, Q), method = "topsoe", unit = "log10")
                  ),
                  sum(((P) * log10(
                          2 * (P) / ((P) + (Q))
                  )) + ((Q) * log10(
                          2 * (Q) / ((P) + (Q))
                  ))))
                  
          })

test_that(
        "distance(method = 'topsoe') computes the correct distance value in case 0 values are stored in the input probability vectors -> leading to 0 * log(0) computations
        .",
        {
                A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                
                topsoe <- function(x, y) {
                        dist <- vector(mode = "numeric", length = 1)
                        dist <- 0
                        
                        for (i in 1:length(x)) {
                                if ((x[i] == 0) & ((y[i]) == 0)) {
                                        dist = dist
                                } else {
                                        dist = dist + (x[i] * log((2 * x[i]) / (x[i] + y[i]))) + (y[i] * log((2 * y[i]) /
                                                                                                                     (x[i] + y[i])))
                                }
                                
                        }
                        
                        return(dist)
                        
                }
                expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "topsoe")), topsoe(A, B))
                
                
        }
)


context("Test implementation of jensen-shannon distance ...")


test_that(
        "distance(method = 'jensen-shannon') computes the correct distance value using unit = log.",
        {
                test_JS_dist <- function(P, Q) {
                        0.5 * ((sum((P) * log((2 * (P)) / ((P) + (Q))
                        )))  +  (sum((Q) * log((2 * (Q)) / ((P) + (Q))
                        ))))
                }
                
                expect_equal(as.vector(
                        philentropy::distance(rbind(P, Q), method = "jensen-shannon")
                ),
                test_JS_dist(P, Q))
                
        }
)


test_that(
        "distance(method = 'jensen-shannon') computes the correct distance value using unit = log2.",
        {
                expect_equal(as.vector(
                        philentropy::distance(rbind(P, Q), method = "jensen-shannon", unit = "log2")
                ),
                0.5 * ((sum((P) * log2((2 * (P)) / ((P) + (Q)))
                ))  +  (sum((Q) * log2((2 * (Q)) / ((P) + (Q)))
                ))))
                
        }
)

test_that(
        "distance(method = 'jensen-shannon') computes the correct distance value using unit = log10.",
        {
                expect_equal(as.vector(
                        philentropy::distance(rbind(P, Q), method = "jensen-shannon", unit = "log10")
                ),
                0.5 * ((sum((P) * log10((2 * (P)) / ((P) + (Q)))
                ))  +  (sum((Q) * log10((2 * (Q)) / ((P) + (Q)))
                ))))
                
        }
)

test_that(
        "distance(method = 'jensen-shannon') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0 * log(0) computation
        .",
        {
                A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                
                js <- function(x, y) {
                        sum1 <- 0
                        sum2 <- 0
                        
                        
                        for (i in 1:length(x)) {
                          PQsum <- x[i] + y[i]
                          
                                if ((x[i] == 0) | ((y[i]) == 0)) {
                                  if (x[i] == 0.0 || PQsum == 0.0) {
                                    sum1  = sum1
                                  } else {
                                    sum1  = sum1 +  (x[i] * log((2.0 * x[i]) / PQsum))
                                  }
                                  if (y[i] == 0.0 || PQsum == 0.0) {
                                    sum2  = sum2
                                  } else {
                                    sum2  = sum2 +  (y[i] * log((2.0 * y[i]) / PQsum))
                                  }
                                } else {
                                        sum1 = sum1 + (x[i] * log((2 * x[i]) / (x[i] + y[i])))
                                        sum2 = sum2 + (y[i] * log((2 * y[i]) / (x[i] + y[i])))
                                }
                                
                        }
                        
                        return(0.5 * (sum1 + sum2))
                        
                }
                expect_equal(as.vector(
                        philentropy::distance(rbind(A, B), method = "jensen-shannon")
                ), js(A, B))
                
                # test that 0 values are treated correctly
                expect_equal(as.vector(
                        philentropy::distance(rbind(c(0.3, 0.3, 0.4), c(0.5, 0.5, 0.0)), method = "jensen-shannon")
                ),
                js(c(0.3, 0.3, 0.4), c(0.5, 0.5, 0.0)))
                
                expect_equal(as.vector(
                        philentropy::distance(rbind(c(0.5, 0.5, 0.0), c(0.5, 0.5, 0.0)), method = "jensen-shannon")
                ),
                js(c(0.5, 0.5, 0.0), c(0.5, 0.5, 0.0)))
                
                expect_equal(as.vector(
                        philentropy::distance(rbind(c(0.3, 0.3, 0.4), c(0.5, 0.5, 0.0)), method = "jensen-shannon")
                ),
                js(c(0.5, 0.5, 0.0), c(0.3, 0.3, 0.4)))
                
                # test correct computation of distance matrix
                A <- c(0.0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- c(0.0, 0.0, 0.25, 0.25, 0.25, 0.25)
                C <- c(0.0, 0.25, 0.0, 0.25, 0.25, 0.25)
                
                distMat <- rbind(A, B, C)
                dist.vals <-
                        distance(distMat, method = "jensen-shannon")
                
                expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                             test_dist_matrix(distMat, FUN = js))
        }
)

context("Test implementation of jensen_difference distance ...")

test_that(
        "distance(method = 'jensen_difference') computes the correct distance value using unit = log.",
        {
                expect_equal(as.vector(
                        philentropy::distance(rbind(P, Q), method = "jensen_difference")
                ),
                sum(((((P) * log((P))) + ((Q) * log((Q)))
                ) / 2) - (((
                        P
                ) + (
                        Q
                )) / 2) * log(((
                        P
                ) + (
                        Q
                )) / 2)))
                
        }
)


test_that(
        "distance(method = 'jensen_difference') computes the correct distance value using unit = log2.",
        {
                expect_equal(as.vector(
                        philentropy::distance(rbind(P, Q), method = "jensen_difference", unit = "log2")
                ),
                sum(((((P) * log2((P))) + ((Q) * log2((Q)))
                ) / 2) - (((
                        P
                ) + (
                        Q
                )) / 2) * log2(((
                        P
                ) + (
                        Q
                )) / 2)))
                
        }
)


test_that(
        "distance(method = 'jensen_difference') computes the correct distance value using unit = log10.",
        {
                expect_equal(as.vector(
                        philentropy::distance(rbind(P, Q), method = "jensen_difference", unit = "log10")
                ),
                sum(((((P) * log10((P))) + ((Q) * log10((Q)))
                ) / 2) - (((
                        P
                ) + (
                        Q
                )) / 2) * log10(((
                        P
                ) + (
                        Q
                )) / 2)))
                
        }
)


test_that(
        "distance(method = 'jensen_difference') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0 * log(0) computation
        .",
        {
                A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                
                js.diff <- function(P, Q) {
                        dist <- vector(mode = "numeric", length = 1)
                        dist <- 0
                        
                        for (i in 1:length(P)) {
                          PQsum <- P[i] + Q[i]
                          if (PQsum == 0.0 || P[i] == 0.0 || Q[i] == 0.0) {
                            if (PQsum == 0.0 && P[i] == 0.0 && Q[i] == 0.0){
                              dist = dist
                            } 
                            if (P[i] == 0.0 && Q[i] > 0.0 && PQsum > 0.0) {
                              dist = dist + ((0.0 + (Q[i] * log(Q[i]))) / 2.0) - ((PQsum / 2.0) * log(PQsum / 2.0)) ;
                            }
                            
                            if (P[i] > 0.0 && Q[i] == 0.0 && PQsum > 0.0) {
                              dist = dist + (((P[i] * log(P[i])) + 0.0 ) / 2.0) - ((PQsum / 2.0) * log(PQsum / 2.0)) ;
                            }
                            
                            if (P[i] == 0.0 && Q[i] == 0.0 && PQsum > 0.0) {
                              dist = dist + 0.0 - ((PQsum / 2.0) * log(PQsum / 2.0)) ;
                            }
                            
                            if (P[i] > 0.0 && Q[i] > 0.0 && PQsum == 0.0) {
                              dist = dist + (((P[i] * log(P[i])) + (Q[i] * log(Q[i]))) / 2.0) - 0.0 ;
                            }
                            
                            if (P[i] > 0.0 && Q[i] == 0.0 && PQsum == 0.0) {
                              dist = dist + (((P[i] * log(P[i])) + 0.0) / 2.0) - 0.0 ;
                            }
                            
                            if (P[i] == 0.0 && Q[i] > 0.0 && PQsum == 0.0) {
                              dist = dist + ((0.0 + (Q[i] * log(Q[i]))) / 2.0) - 0.0 ;
                            }
                            
                          } else {
                            dist = dist + (((P[i] * log(P[i])) + (Q[i] * log(Q[i]))) / 2.0) - ((PQsum / 2.0) * log(PQsum / 2.0)) ;
                          }   
                                
                        }
                        
                        return(dist)
                        
                }
                expect_equal(as.vector(
                        philentropy::distance(rbind(A, B), method = "jensen_difference")
                ), js.diff(A, B))
                
        }
)


context("Test implementation of taneja distance ...")


test_that("distance(method = 'taneja') computes the correct distance value using unit = log.",
          {
                  test_taneja_dist <- function(P, Q) {
                          sum(((P + Q) / 2) * log((P + Q) / (2 * sqrt(
                                  P * Q
                          ))))
                  }
                  
                  expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "taneja")),
                               test_taneja_dist(P, Q))
                  
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <- distance(distMat, method = "taneja")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_taneja_dist))
                  
          })

test_that("distance(method = 'taneja') computes the correct distance value using unit = log2.",
          {
                  expect_equal(as.vector(
                          philentropy::distance(rbind(P, Q), method = "taneja", unit = "log2")
                  ),
                  sum(((P + Q) / 2) * log2((P + Q) / (
                          2 * sqrt(P * Q)
                  ))))
                  
          })

test_that("distance(method = 'taneja') computes the correct distance value using unit = log10.",
          {
                  expect_equal(as.vector(
                          philentropy::distance(rbind(P, Q), method = "taneja", unit = "log10")
                  ),
                  sum(((P + Q) / 2) * log10((P + Q) / (
                          2 * sqrt(P * Q)
                  ))))
                  
          })

test_that(
        "distance(method = 'taneja') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0 * log(0) computation
        .",
        {
                A <- c(0, 0.25, 0.25, 0, 0.25, 0.25)
                B <- c(0, 0, 0.25, 0.25, 0.25, 0.25)
                
                taneja <- function(x, y) {
                        dist <- vector(mode = "numeric", length = 1)
                        dist <- 0
                        
                        for (i in 1:length(x)) {
                                if ((x[i] == 0) & ((y[i]) == 0)) {
                                        dist = dist
                                } else {
                                        denominator <- (2 * sqrt(x[i] * y[i]))
                                        
                                        if (denominator == 0) {
                                                dist = dist + (((
                                                        x[i] + y[i]
                                                ) / 2) * log((
                                                        x[i] + y[i]
                                                ) / 0.00001))
                                        } else {
                                                dist = dist + (((
                                                        x[i] + y[i]
                                                ) / 2) * log((x[i] + y[i]) / denominator
                                                ))
                                        }
                                }
                        }
                        
                        return(dist)
                        
                }
                expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "taneja")), taneja(A, B))
                
        }
)


context("Test implementation of kumar-johnson distance ...")

test_that("distance(method = 'kumar-johnson') computes the correct distance value.",
          {
                  test_kumar_dist <- function(P, Q) {
                          sum(((P ^ 2 - Q ^ 2) ^ 2 / (2 * (
                                  P * Q
                          ) ^ 1.5)))
                  }
                  
                  expect_equal(as.vector(
                          philentropy::distance(rbind(P, Q), method = "kumar-johnson")
                  ),
                  test_kumar_dist(P, Q))
                  
                  # test correct computation of distance matrix
                  distMat <-
                          rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
                  dist.vals <-
                          distance(distMat, method = "kumar-johnson")
                  
                  expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                               test_dist_matrix(distMat, FUN = test_kumar_dist))
          })


context("Test implementation of avg distance ...")

test_that("distance(method = 'avg') computes the correct distance value.", {
        test_avg_dist <- function(P, Q) {
                (sum(abs(P - Q)) + max(abs(P - Q))) / 2
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "avg")),
                     test_avg_dist(P, Q))
        
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2, 5), rep(0.1, 5), c(5, 1, 7, 9, 5))
        dist.vals <- distance(distMat, method = "avg")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_avg_dist))
        
})
