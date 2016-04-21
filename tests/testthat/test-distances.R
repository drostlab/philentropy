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



P <- 1:10/sum(1:10)
Q <- 20:29/sum(20:29)
V <- -10:10
W <- -20:0

# function to test distance matrix functionality
# for different distance measures
test_dist_matrix <- function(x, FUN){
        
        dist.fun <- match.fun(FUN)
        res.dist.matrix <- matrix(NA_real_,nrow(x),nrow(x))
        
        for(i in 1:nrow(x)){
                for(j in 1:nrow(x)){
                        res.dist.matrix[i,j] <- dist.fun(x[i, ],x[j, ])
                }
        }
        return(res.dist.matrix[lower.tri(res.dist.matrix, diag = FALSE)])
}

test_that("'euclidien' (or any other not implemented string) is caught when wrong input string for method is entered", {
        
        distMat <- rbind(rep(0.2,5),rep(0.1,10))
        expect_error(distance(distMat, method = "euclidien"),"Method 'euclidien' is not implemented in this function. Please consult getDistMethods().")
} )

test_that("Only numeric values are passed to distance()", {

        distMat <- rbind(rep("A",10),rep("B",10))
        expect_error(distance(distMat, method = "euclidean"), paste0("Your input ",class(distMat)," stores non-numeric values. Non numeric values cannot be used to compute distances.."))
})


test_that("Only choose from units: log, log2, or log10", {
        
        distMat <- rbind(rep(0.2,5),rep(0.1,10))
        expect_error(distance(distMat, method = "euclidean", unit = "log5"), 
                     "You can only choose units: log, log2, or log10.")
})

test_that("distance(method = 'euclidean') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "euclidean")),
                     sqrt(sum(abs((P) - (Q))^2)))
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "euclidean")), 
                     as.vector(stats::dist(base::rbind(P,Q), method = "euclidean")))
        
        expect_equal(as.vector(philentropy::distance(rbind(V, W), method = "euclidean")), 
                     as.vector(stats::dist(base::rbind(V,W), method = "euclidean")))
        
        expect_equal(as.vector(philentropy::distance(rbind(V, W), method = "euclidean")),
                     sqrt(sum(abs((V) - (W))^2)))
        
        
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "euclidean")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     as.vector(dist(distMat)))
        
        #expect_error(philentropy::distance(1:10, 20:29, method = "euclidean"))
})


test_that("distance(method = 'manhattan') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "manhattan")),
                     sum(abs((P) - (Q))))
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "manhattan")),
                     as.vector(stats::dist(base::rbind(P,Q), method = "manhattan")))
        
        expect_equal(as.vector(philentropy::distance(rbind(V, W), method = "manhattan")),
                     sum(abs((V) - (W))))
        
        expect_equal(as.vector(philentropy::distance(rbind(V, W), method = "manhattan")), 
                     as.vector(stats::dist(base::rbind(V,W), method = "manhattan")))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "manhattan")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     as.vector(dist(distMat, method = "manhattan")))
})


test_that("distance(method = 'minkowski') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "minkowski", p = 4)),
                     (sum(abs((P) - (Q))^4))^0.25)
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "minkowski", p = 4)), 
                     as.vector(stats::dist(base::rbind(P,Q), method = "minkowski", p = 4)))
        
        expect_error(as.vector(philentropy::distance(rbind(P, Q), method = "minkowski")),
                     "Please specify p for the Minkowski distance!")
        
        expect_equal(as.vector(philentropy::distance(rbind(V, W), method = "minkowski", p = 4)), (sum(abs((V) - (W))^4))^0.25)
        
        expect_equal(as.vector(philentropy::distance(rbind(V, W), method = "minkowski", p = 4)),
                     as.vector(stats::dist(base::rbind(V,W), method = "minkowski", p = 4)))
        
        expect_error(as.vector(philentropy::distance(rbind(V, W), method = "minkowski")),
                     "Please specify p for the Minkowski distance!")
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "minkowski", p = 4)
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     as.vector(dist(distMat, method = "minkowski", p = 4)))
        
})


test_that("distance(method = 'chebyshev') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "chebyshev")),
                     max(abs((P) - (Q))))
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "chebyshev")),
                     as.vector(stats::dist(base::rbind(P,Q), method = "maximum")))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "chebyshev")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     as.vector(dist(distMat, method = "maximum")))
        
})


test_that("distance(method = 'sorensen') computes the correct distance value.", {
        
        test_sorensen_dist <- function(P,Q){
                sum(abs((P) - (Q))) / sum((P) + (Q))
        }
        
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "sorensen")),
                     test_sorensen_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "sorensen")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_sorensen_dist))
        
})


test_that("distance(method = 'gower') computes the correct distance value.", {
        
        test_gower_dist <- function(P,Q){
                return((1/length(P)) * sum(abs((P) - (Q))))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "gower")), 
                     test_gower_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "gower")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_gower_dist))
        
})


test_that("distance(method = 'soergel') computes the correct distance value.", {
        
        test_soergel_dist <- function(P,Q){
                sum(abs((P) - (Q))) / sum(apply(rbind(P, Q),2,max))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "soergel")), 
                     test_soergel_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "soergel")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_soergel_dist))
        
})


test_that("distance(method = 'kulczynski_d') computes the correct distance value.", {
        
        test_kd_dist <- function(P,Q){
                sum(abs((P) - (Q))) / sum(apply(rbind(P, Q),2,min))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "kulczynski_d")),
                     test_kd_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "kulczynski_d")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_kd_dist))
        
})


test_that("distance(method = 'canberra') computes the correct distance value.", {
        
        test_canberra_dist <- function(P,Q){
                sum( abs((P) - (Q)) / ((P) + (Q)))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "canberra")),
                     test_canberra_dist(P,Q))
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "canberra")),
                     as.vector(stats::dist(base::rbind(P,Q), method = "canberra")))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "canberra")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_canberra_dist))
        
})


test_that("distance(method = 'canberra') computes the correct distance value when P_i and Q_i are 0 -> 0/0 is then replaced by 0.", {
        
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        
        canb <- function(x,y){
                
                dist <- vector(mode = "numeric", length = 1)
                dist <- 0
                
                for(i in 1:length(x)){
                        
                        if((abs(x[i] - y[i]) == 0) & ((x[i] + y[i]) == 0)){
                                dist = dist
                        } else {
                                dist = dist + (abs(x[i] - y[i]) / ((x[i]) + (y[i])))
                        }
                        
                }
                
                return(dist)
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "canberra")), canb(A,B))
        
        
})

test_that("distance(method = 'lorentzian') computes the correct distance value using unit = log.", {
        
        test_lorentzian_dist <- function(P,Q){
                sum( log(1 + abs((P) - (Q))))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "lorentzian")),
                     test_lorentzian_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "lorentzian")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_lorentzian_dist))
        
})

test_that("distance(method = 'lorentzian') computes the correct distance value using unit = log2.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "lorentzian",
                                                     unit = "log2")), sum( log2(1 + abs((P) - (Q)))))
        
})

test_that("distance(method = 'lorentzian') computes the correct distance value using unit = log10.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "lorentzian",
                                                     unit = "log10")), sum( log10(1 + abs((P) - (Q)))))
        
        
})


test_that("distance(method = 'intersection') computes the correct distance value.", {
        
        test_intersection_dist <- function(P,Q){
                sum(apply(base::rbind(P,Q),2,min))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "intersection")),
                     test_intersection_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "intersection")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_intersection_dist))
        
})


test_that("distance(method = 'non-intersection') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "non-intersection")),
                     1 - sum(apply(base::rbind(P,Q),2,min)))
        
})


test_that("distance(method = 'wavehedges') computes the correct distance value.", {
        
        test_wh_dist <- function(P,Q){
                sum(abs(P - Q) / apply(base::rbind(P,Q),2,max))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "wavehedges")),
                     test_wh_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "wavehedges")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_wh_dist))
})



test_that("distance(method = 'wavehedges') computes the correct distance value in case the input probability vectors store 0 values at the same position causing 0/0 computations.", {
        
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        
        wh <- function(x,y){
                
                dist <- vector(mode = "numeric", length = 1)
                dist <- 0
                
                for(i in 1:length(x)){
                        
                        if((abs(x[i] - y[i]) == 0) & ((max(c(x[i],y[i]))) == 0)){
                                dist = dist
                        } else {
                                dist = dist + (abs(x[i] - y[i]) / max(c(x[i],y[i])))
                        }
                }
                return(dist)
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "wavehedges")), wh(A,B))
        
})


test_that("distance(method = 'czekanowski') computes the correct distance value.", {
        
        test_czekanowski_dist <- function(P,Q){
                sum(abs(P - Q)) / sum(P + Q)
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "czekanowski")),
                     test_czekanowski_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "czekanowski")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_czekanowski_dist))
        
})


test_that("distance(method = 'motyka') computes the correct distance value.", {
        
        test_motyka_dist <- function(P,Q){
                sum(apply(base::rbind(P,Q),2,max)) / sum(P + Q)
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "motyka")), 
                     test_motyka_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "motyka")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_motyka_dist))
})


test_that("distance(method = 'kulczynski_s') computes the correct distance value.", {
        
        test_ks_dist <- function(P,Q){
                1 / (sum(abs((P) - (Q))) / sum(apply(rbind(P, Q),2,min)))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "kulczynski_s")),
                     test_ks_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "kulczynski_s")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_ks_dist))
        
})


test_that("distance(method = 'tanimoto') computes the correct distance value.", {
        
        test_tanimoto_dist <- function(P,Q){
                sum(abs((P) - (Q))) / sum(apply(rbind(P, Q),2,max))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "tanimoto")),
                     test_tanimoto_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "tanimoto")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_tanimoto_dist))
})


test_that("distance(method = 'ruzicka') computes the correct distance value.", {
        
        test_ruzicka_dist <- function(P,Q){
                1 - (sum(abs((P) - (Q))) / sum(apply(rbind(P, Q),2,max)))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "ruzicka")),
                     test_ruzicka_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "ruzicka")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_ruzicka_dist))
        
})


test_that("distance(method = 'inner_product') computes the correct distance value.", {
        
        test_innerproduct_dist <- function(P,Q){
                sum ( P*Q )
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "inner_product")),
                     test_innerproduct_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "inner_product")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_innerproduct_dist))
})


test_that("distance(method = 'harmonic_mean') computes the correct distance value.", {
        
        test_harmonic_mean_dist <- function(P,Q){
                2 * sum ( (P * Q) / (P + Q) )
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "harmonic_mean")),
                     test_harmonic_mean_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "harmonic_mean")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_harmonic_mean_dist))
        
})


test_that("distance(method = 'harmonic_mean') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0/0 computation

.", {
        
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        
        hm <- function(x,y){
                
                dist <- vector(mode = "numeric", length = 1)
                dist <- 0
                
                for(i in 1:length(x)){
                        
                        if(((x[i] * y[i]) == 0) & ((x[i] + y[i]) == 0)){
                                dist = dist
                        } else {
                                dist = dist + ((x[i] * y[i]) / (x[i] + y[i]))
                        }
                        
                }
                
                return(2 * dist)
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "harmonic_mean")), hm(A,B))        
})


test_that("distance(method = 'cosine') computes the correct distance value.", {
        
        test_cosine_dist <- function(P,Q){
                ((sum ( (P) * (Q) )) / (sqrt(sum((P)^2)) * sqrt(sum((Q)^2))))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "cosine")),
                     test_cosine_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "cosine")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_cosine_dist))
        
})


test_that("distance(method = 'hassebrook') computes the correct distance value.", {
        
        test_hassebrook_dist <- function(P,Q){
                ((sum ( (P) * (Q) )) / (sum((P)^2) + sum((Q)^2) - ((sum ( (P) * (Q) )))))
        }
        
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "hassebrook")),
                     test_hassebrook_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "hassebrook")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_hassebrook_dist))
        
})



test_that("distance(method = 'jaccard') computes the correct distance value.", {
        
        test_jaccard_dist <- function(P,Q){
                1 - ((sum ( (P) * (Q) )) / (sum((P)^2) + sum((Q)^2) - ((sum ( (P) * (Q) )))))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "jaccard")),
                     test_jaccard_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "jaccard")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_jaccard_dist))
})


test_that("distance(method = 'dice') computes the correct distance value.", {
        
        test_dice_dist <- function(P,Q){
                1 - (2 * (sum ( (P) * (Q) )) / (sum((P)^2) + sum((Q)^2) ))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "dice")),
                     test_dice_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "dice")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_dice_dist))
        
})


test_that("distance(method = 'fidelity') computes the correct distance value.", {
        
        test_fidelity_dist <- function(P,Q){
                sum(sqrt(P * Q))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "fidelity")),
                     test_fidelity_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "fidelity")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_fidelity_dist))
        
})


test_that("distance(method = 'bhattacharyya') computes the correct distance value using unit = log.", {
        
        test_b_dist <- function(P,Q){
                -log(sum(sqrt(P * Q)))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "bhattacharyya")),
                     test_b_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "bhattacharyya")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_b_dist))
        
})

test_that("distance(method = 'bhattacharyya') computes the correct distance value using unit = log2.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "bhattacharyya", unit = "log2")), -log2(sum(sqrt(P * Q))))
        
})

test_that("distance(method = 'bhattacharyya') computes the correct distance value using unit = log10.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "bhattacharyya", unit = "log10")), -log10(sum(sqrt(P * Q))))
        
})


test_that("distance(method = 'hellinger') computes the correct distance value.", {
        
        test_hellinger_dist <- function(P,Q){
                2L * sqrt(1L - sum(sqrt(P * Q)))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "hellinger")),
                     test_hellinger_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5)/sum(c(5,1,7,9,5)))
        dist.vals <- distance(distMat, method = "hellinger")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_hellinger_dist))
        
})



test_that("distance(method = 'matusita') computes the correct distance value.", {
        
        test_matusita_dist <- function(P,Q){
                sqrt(sum((sqrt(P) - sqrt(Q))^2))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "matusita")),
                     test_matusita_dist(P,Q))
        
        
        # test correct computation of distance matrix
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        C <- c(0,0.25,0,0.25,0.25,0.25)
        distMat <- rbind(A,B,C)
        dist.vals <- distance(distMat, method = "matusita")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_matusita_dist))
})


test_that("distance(method = 'squared_chord') computes the correct distance value.", {
        
        test_sqchord_dist <- function(P,Q){
                sum((sqrt(P) - sqrt(Q))^2)
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "squared_chord")),
                     test_sqchord_dist(P,Q))
        
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "squared_chord")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_sqchord_dist))
})


test_that("distance(method = 'squared_euclidean') computes the correct distance value.", {
        
        test_sqeuclidean_dist <- function(P,Q){
                sum(((P) - (Q))^2)
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "squared_euclidean")), 
                     test_sqeuclidean_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "squared_euclidean")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_sqeuclidean_dist))
        
})


test_that("distance(method = 'pearson') computes the correct distance value.", {
        
        test_pearson_dist <- function(P,Q){
                sum((P - Q)^2 / Q)
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "pearson")),
                     test_pearson_dist(P,Q))
        
        # test correct computation of distance matrix
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        C <- c(0,0.25,0,0.25,0.25,0.25)
        distMat <- rbind(A,B,C)
        dist.vals <- distance(distMat, method = "pearson")
        
#         expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
#                      test_dist_matrix(distMat, FUN = test_pearson_dist))
        
})


test_that("distance(method = 'neyman') computes the correct distance value.", {
        
        test_neyman_dist <- function(P,Q){
                sum(((P - Q)^2) / P)
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "neyman")),
                     test_neyman_dist(P,Q))
        
        # test correct computation of distance matrix
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        C <- c(0,0.25,0,0.25,0.25,0.25)
        distMat <- rbind(A,B,C)
        dist.vals <- distance(distMat, method = "neyman")
        
#         expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
#                      test_dist_matrix(distMat, FUN = test_neyman_dist))
})


test_that("distance(method = 'squared_chi') computes the correct distance value.", {
        
        test_sqchi_dist <- function(P,Q){
                sum(((P) - (Q))^2 / ((P) + (Q)))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "squared_chi")),
                     test_sqchi_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "squared_chi")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_sqchi_dist))
})



test_that("distance(method = 'squared_chi') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0/0 computation
.", {
        
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        
        sqchisq <- function(x,y){
                
                dist <- vector(mode = "numeric", length = 1)
                dist <- 0
                
                for(i in 1:length(x)){
                        
                        if(((x[i] - y[i])^2 == 0) & ((x[i] + y[i]) == 0)){
                                dist = dist
                        } else {
                                dist = dist + ((x[i] - y[i])^2 / (x[i] + y[i]))
                        }
                        
                }
                
                return(dist)
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "squared_chi")), sqchisq(A,B))
        
})



test_that("distance(method = 'prob_symm') computes the correct distance value.", {
        
        test_probsymm_dist <- function(P,Q){
                2 * sum(((P) - (Q))^2 / ((P) + (Q)))
        }  
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "prob_symm")),
                     test_probsymm_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "prob_symm")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_probsymm_dist))
})


test_that("distance(method = 'prob_symm') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0/0 computation
.", {
        
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        
        probsymmchisq <- function(x,y){
                
                dist <- vector(mode = "numeric", length = 1)
                dist <- 0
                
                for(i in 1:length(x)){
                        
                        if(((x[i] - y[i])^2 == 0) & ((x[i] + y[i]) == 0)){
                                dist = dist
                        } else {
                                dist = dist + ((x[i] - y[i])^2 / (x[i] + y[i]))
                        }
                        
                }
                
                return(2 * dist)
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "prob_symm")), probsymmchisq(A,B))
        
})

test_that("distance(method = 'divergence') computes the correct distance value.", {
        
        test_divergence_dist <- function(P,Q){
                2 * sum(((P) - (Q))^2 / ((P) + (Q))^2)
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "divergence")),
                     test_divergence_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "divergence")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_divergence_dist))
})


test_that("distance(method = 'divergence') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0/0 computation
.", {
        
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        
        div <- function(x,y){
                
                dist <- vector(mode = "numeric", length = 1)
                dist <- 0
                
                for(i in 1:length(x)){
                        
                        if(((x[i] - y[i])^2 == 0) & ((x[i] + y[i])^2 == 0)){
                                dist = dist
                        } else {
                                dist = dist + ((x[i] - y[i])^2 / (x[i] + y[i])^2)
                        }
                        
                }
                
                return(2 * dist)
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "divergence")), div(A,B))
})

test_that("distance(method = 'clark') computes the correct distance value.", {
        
        test_clark_dist <- function(P,Q){
                sqrt(sum((abs((P) - (Q)) / ((P) + (Q)))^2))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "clark")),
                     test_clark_dist(P,Q))
        
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "clark")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_clark_dist))
})


test_that("distance(method = 'clark') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0/0 computation
.", {
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        
        clark <- function(x,y){
                
                dist <- vector(mode = "numeric", length = 1)
                dist <- 0
                
                for(i in 1:length(x)){
                        
                        if((abs(x[i] - y[i]) == 0) & ((x[i] + y[i]) == 0)){
                                dist = dist
                        } else {
                                dist = dist + (abs(x[i] - y[i]) / (x[i] + y[i]))^2
                        }
                        
                }
                
                return(sqrt(dist))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "clark")), clark(A,B))
        
})


test_that("distance(method = 'additive_symm') computes the correct distance value.", {
        
        test_addsymm_dist <- function(P,Q){
                sum((((P) - (Q))^2 * ((P) + (Q))) / ((P) * (Q)))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "additive_symm")),
                     test_addsymm_dist(P,Q))
        
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "additive_symm")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_addsymm_dist))
})


test_that("distance(method = 'additive_symm') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0/0 computation
.", {
        
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        
        add <- function(x,y){
                
                dist <- vector(mode = "numeric", length = 1)
                dist <- 0
                
                for(i in 1:length(x)){
                        
                        if(((x[i] + y[i]) == 0) & ((x[i] * y[i]) == 0)){
                                dist = dist
                        } else {
                                dist = dist + ((x[i] - y[i])^2 * ((x[i] + y[i]) / (x[i] * y[i])))
                        }
                        
                }
                
                return(dist)
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "additive_symm")), add(A,B))
})

test_that("distance(method = 'kullback-leibler') computes the correct distance value using unit = log.", {
        
        test_KL_dist <- function(P,Q){
                sum((P) * log((P) / (Q)))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "kullback-leibler")),
                     test_KL_dist(P,Q))
        
        # test correct computation of distance matrix
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        C <- c(0,0.25,0,0.25,0.25,0.25)
        distMat <- rbind(A,B,C)
        dist.vals <- distance(distMat, method = "kullback-leibler")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_KL_dist))
})


test_that("distance(method = 'kullback-leibler') computes the correct distance value using unit = log2.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "kullback-leibler", unit = "log2")), sum((P) * log2((P) / (Q))))
        
})

test_that("distance(method = 'kullback-leibler') computes the correct distance valueusing unit = log10.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "kullback-leibler", unit = "log10")), sum((P) * log10((P) / (Q))))
        
})

test_that("distance(method = 'kullback-leibler') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0 * log(0) computation
.", {

        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        
        kl <- function(x,y){
                
                dist <- vector(mode = "numeric", length = 1)
                dist <- 0
                
                for(i in 1:length(x)){
                        
                        if((x[i] == 0) & ((y[i]) == 0)){
                                dist = dist
                        } else {
                                dist = dist + (x[i] * log(x[i]/y[i]))
                        }
                        
                }
                
                return(dist)
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "kullback-leibler")), kl(A,B))
        
})




test_that("distance(method = 'jeffreys') computes the correct distance value using unit = log.", {
        
        test_jeffreys_dist <- function(P,Q){
                sum(((P) - (Q)) * log((P) / (Q)))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "jeffreys")),
                     test_jeffreys_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "jeffreys")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_jeffreys_dist))
        
})


test_that("distance(method = 'jeffreys') computes the correct distance value using unit = log2.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "jeffreys", unit = "log2")), sum(((P) - (Q)) * log2((P) / (Q))))
        
})

test_that("distance(method = 'jeffreys') computes the correct distance value using unit = log10.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "jeffreys", unit = "log10")), sum(((P) - (Q)) * log10((P) / (Q))))
        
})

test_that("distance(method = 'jeffreys') computes the correct distance value in case 0 values are stored in the input probability vectors -> leading to log(0) computations.", {
        
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- rep(1,6) / 6
        
        jeff <- function(x,y){
                
                if(any((x/y) == 0)){
                        xy.ratio <- x/y
                        xy.ratio[xy.ratio == 0] <- 0.00001
                        sum((x - y) * log(xy.ratio))
                } else {
                        
                        sum((x - y) * log(x/y))
                }
                
        }
        expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "jeffreys")), jeff(A,B))
})



test_that("distance(method = 'k_divergence') computes the correct distance value using unit = log.", {
        
        test_kdivergence_dist <- function(P,Q){
                sum((P) * log(2 * (P) / ((P) + (Q))))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "k_divergence")),
                     test_kdivergence_dist(P,Q))
        
        # test correct computation of distance matrix
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        C <- c(0,0.25,0,0.25,0.25,0.25)
                
        distMat <- rbind(A,B,C)
        dist.vals <- distance(distMat, method = "k_divergence")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_kdivergence_dist))
})

test_that("distance(method = 'k_divergence') computes the correct distance value using unit = log2.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "k_divergence", unit = "log2")), sum((P) * log2(2 * (P) / ((P) + (Q)))))
        
})

test_that("distance(method = 'k_divergence') computes the correct distance value using unit = log10.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "k_divergence", unit = "log10")), sum((P) * log10(2 * (P) / ((P) + (Q)))))
        
})

test_that("distance(method = 'k_divergence') computes the correct distance value in case 0 values are stored in the input probability vectors -> leading to 0 * log(0) computations.", {
        
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        
        kdiv <- function(x,y){
                
                dist <- vector(mode = "numeric", length = 1)
                dist <- 0
                
                for(i in 1:length(x)){
                        
                        if((x[i] == 0) & ((y[i]) == 0)){
                                dist = dist
                        } else {
                                dist = dist + (x[i] * log((2 * x[i])/(x[i]+y[i])))
                        }
                        
                }
                
                return(dist)
                
        }
        expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "k_divergence")), kdiv(A,B))
})


test_that("distance(method = 'topsoe') computes the correct distance value using unit = log.", {
        
        test_topsoe_dist <- function(P,Q){
                sum(((P) * log(2 * (P) / ((P) + (Q)))) + ((Q) * log(2 * (Q) / ((P) + (Q)))))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "topsoe")),
                     test_topsoe_dist(P,Q))
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "topsoe")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_topsoe_dist))
        
})


test_that("distance(method = 'topsoe') computes the correct distance value using unit = log2.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "topsoe", unit = "log2")),
                     sum(((P) * log2(2 * (P) / ((P) + (Q)))) + ((Q) * log2(2 * (Q) / ((P) + (Q))))))
        
})


test_that("distance(method = 'topsoe') computes the correct distance value using unit = log10.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "topsoe", unit = "log10")),
                     sum(((P) * log10(2 * (P) / ((P) + (Q)))) + ((Q) * log10(2 * (Q) / ((P) + (Q))))))
        
})

test_that("distance(method = 'topsoe') computes the correct distance value in case 0 values are stored in the input probability vectors -> leading to 0 * log(0) computations
.", {
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        
        topsoe <- function(x,y){
                
                dist <- vector(mode = "numeric", length = 1)
                dist <- 0
                
                for(i in 1:length(x)){
                        
                        if((x[i] == 0) & ((y[i]) == 0)){
                                dist = dist
                        } else {
                                dist = dist + (x[i] * log((2 * x[i])/(x[i]+y[i]))) + (y[i] * log((2 * y[i])/(x[i]+y[i])))
                        }
                        
                }
                
                return(dist)
                
        }
        expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "topsoe")), topsoe(A,B))
        
        
})

test_that("distance(method = 'jensen-shannon') computes the correct distance value using unit = log.", {
        
        test_JS_dist <- function(P,Q){
                0.5 * ((sum((P) * log((2 * (P)) / ((P) + (Q)))))  +  (sum((Q) * log((2 * (Q)) / ((P) + (Q))))))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "jensen-shannon")),
                     test_JS_dist(P,Q))
        
        # test correct computation of distance matrix
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        C <- c(0,0.25,0,0.25,0.25,0.25)
        
        distMat <- rbind(A,B,C)
        dist.vals <- distance(distMat, method = "jensen-shannon")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_JS_dist))
})


test_that("distance(method = 'jensen-shannon') computes the correct distance value using unit = log2.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "jensen-shannon", unit = "log2")),
                     0.5 * ((sum((P) * log2((2 * (P)) / ((P) + (Q)))))  +  (sum((Q) * log2((2 * (Q)) / ((P) + (Q)))))))
        
})

test_that("distance(method = 'jensen-shannon') computes the correct distance value using unit = log10.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "jensen-shannon", unit = "log10")),
                     0.5 * ((sum((P) * log10((2 * (P)) / ((P) + (Q)))))  +  (sum((Q) * log10((2 * (Q)) / ((P) + (Q)))))))
        
})

test_that("distance(method = 'jensen-shannon') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0 * log(0) computation
.", {
        
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        
        js <- function(x,y){
                
                sum1 <- 0
                sum2 <- 0
                
                for(i in 1:length(x)){
                        
                        if((x[i] == 0) & ((y[i]) == 0)){
                                sum1 = sum1
                                sum2 = sum2
                        } else {
                                sum1 = sum1 + (x[i] * log((2 * x[i])/(x[i]+y[i])))
                                sum2 = sum2 + (y[i] * log((2 * y[i])/(x[i]+y[i])))
                        }
                        
                }
                
                return(0.5 * (sum1 + sum2))
                
        }
        expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "jensen-shannon")), js(A,B))
})


test_that("distance(method = 'jensen_difference') computes the correct distance value using unit = log.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "jensen_difference")),
                     sum(((((P) * log((P))) + ((Q) * log((Q)))) / 2 ) - (((P) + (Q)) / 2) * log(((P) + (Q)) / 2)))
        
})


test_that("distance(method = 'jensen_difference') computes the correct distance value using unit = log2.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "jensen_difference", unit = "log2")),
                     sum(((((P) * log2((P))) + ((Q) * log2((Q)))) / 2 ) - (((P) + (Q)) / 2) * log2(((P) + (Q)) / 2)))
        
})


test_that("distance(method = 'jensen_difference') computes the correct distance value using unit = log10.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "jensen_difference", unit = "log10")),
                     sum(((((P) * log10((P))) + ((Q) * log10((Q)))) / 2 ) - (((P) + (Q)) / 2) * log10(((P) + (Q)) / 2)))
        
})


test_that("distance(method = 'jensen_difference') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0 * log(0) computation
.", {
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        
        js.diff <- function(x,y){
                
                dist <- vector(mode = "numeric", length = 1)
                dist <- 0
                
                for(i in 1:length(x)){
                        
                        if((x[i] == 0) & ((y[i]) == 0)){
                                dist = dist
                        } else {
                                dist = dist + ((((x[i] * log(x[i])) + (y[i] * log(y[i]))) / 2) - (((x[i] + y[i])/2) * log(((x[i] + y[i])/2))))
                                
                        }
                        
                }
                
                return(dist)
                
        }
        expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "jensen_difference")), js.diff(A,B))
       
})

test_that("distance(method = 'taneja') computes the correct distance value using unit = log.", {
        
        test_taneja_dist <- function(P,Q){
                sum(((P + Q) / 2) * log((P+Q) / (2 * sqrt(P*Q))))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "taneja")),
                     test_taneja_dist(P,Q) )
        
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "taneja")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_taneja_dist))
        
})

test_that("distance(method = 'taneja') computes the correct distance value using unit = log2.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "taneja", unit = "log2")),
                     sum(((P + Q) / 2) * log2((P+Q) / (2 * sqrt(P*Q)))) )
        
})

test_that("distance(method = 'taneja') computes the correct distance value using unit = log10.", {
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "taneja", unit = "log10")),
                     sum(((P + Q) / 2) * log10((P+Q) / (2 * sqrt(P*Q)))) )
        
})

test_that("distance(method = 'taneja') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0 * log(0) computation
.", {
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        
        taneja <- function(x,y){
                
                dist <- vector(mode = "numeric", length = 1)
                dist <- 0
                
                for(i in 1:length(x)){
                        
                        if((x[i] == 0) & ((y[i]) == 0)){
                                dist = dist
                        } else {
                                
                                denominator <- (2 * sqrt(x[i] * y[i]))
                                
                                if(denominator == 0){
                                        dist = dist + (((x[i] + y[i])/2) * log((x[i] + y[i]) / 0.00001))
                                } else {
                                        dist = dist + (((x[i] + y[i])/2) * log((x[i] + y[i]) / denominator))
                                }
                        }
                }
                
                return(dist)
                
        }
        expect_equal(as.vector(philentropy::distance(rbind(A, B), method = "taneja")), taneja(A,B))
        
})

test_that("distance(method = 'kumar-johnson') computes the correct distance value.", {
        
        test_kumar_dist <- function(P,Q){
                sum(((P^2 - Q^2)^2 / (2 * (P*Q)^1.5)))
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "kumar-johnson")),
                     test_kumar_dist(P,Q) )
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "kumar-johnson")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_kumar_dist))
})



test_that("distance(method = 'avg') computes the correct distance value.", {
        
        test_avg_dist <- function(P,Q){
                (sum( abs(P-Q) ) + max( abs(P-Q))) / 2
        }
        
        expect_equal(as.vector(philentropy::distance(rbind(P, Q), method = "avg")),
                     test_avg_dist(P,Q) )
        
        
        # test correct computation of distance matrix
        distMat <- rbind(rep(0.2,5),rep(0.1,5), c(5,1,7,9,5))
        dist.vals <- distance(distMat, method = "avg")
        
        expect_equal(dist.vals[lower.tri(dist.vals, diag = FALSE)],
                     test_dist_matrix(distMat, FUN = test_avg_dist))
        
})

