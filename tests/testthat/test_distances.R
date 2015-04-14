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


P <- 1:10/sum(1:10)
Q <- 20:29/sum(20:29)

context("Test implementation of distance measures...")

test_that("distance(method = 'euclidean') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "euclidean")), sqrt(sum(abs((P) - (Q))^2)))
        expect_equal(as.vector(philentropy::distance(P, Q, method = "euclidean")), as.vector(stats::dist(base::rbind(P,Q), method = "euclidean")))
        #expect_error(philentropy::distance(1:10, 20:29, method = "euclidean"))
})


test_that("distance(method = 'manhattan') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "manhattan")), sum(abs((P) - (Q))))
        expect_equal(as.vector(philentropy::distance(P, Q, method = "manhattan")), as.vector(stats::dist(base::rbind(P,Q), method = "manhattan")))
})


test_that("distance(method = 'minkowski') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "minkowski", p = 4)), (sum(abs((P) - (Q))^4))^0.25)
        expect_equal(as.vector(philentropy::distance(P, Q, method = "minkowski", p = 4)), as.vector(stats::dist(base::rbind(P,Q), method = "minkowski", p = 4)))
        
})


test_that("distance(method = 'chebyshev') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "chebyshev")), max(abs((P) - (Q))))
        expect_equal(as.vector(philentropy::distance(P, Q, method = "chebyshev")), as.vector(stats::dist(base::rbind(P,Q), method = "maximum")))
        
})


test_that("distance(method = 'sorensen') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "sorensen")), sum(abs((P) - (Q))) / sum((P) + (Q)))
        
})


test_that("distance(method = 'gower') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "gower")), (1/length(1:10)) * sum(abs((P) - (Q))))
        
})


test_that("distance(method = 'soergel') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "soergel")), sum(abs((P) - (Q))) / sum(apply(rbind(P, Q),2,max)))
        
})


test_that("distance(method = 'kulczynski_d') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "kulczynski_d")), sum(abs((P) - (Q))) / sum(apply(rbind(P, Q),2,min)))
        
})


test_that("distance(method = 'canberra') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "canberra")), sum( abs((P) - (Q)) / ((P) + (Q))))
        expect_equal(as.vector(philentropy::distance(P, Q, method = "canberra")), as.vector(stats::dist(base::rbind(P,Q), method = "canberra")))
        
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
        
        expect_equal(as.vector(philentropy::distance(A, B, method = "canberra")), canb(A,B))
        
        
})

test_that("distance(method = 'lorentzian') computes the correct distance value using unit = log.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "lorentzian")), sum( log(1 + abs((P) - (Q)))))
        
})

test_that("distance(method = 'lorentzian') computes the correct distance value using unit = log2.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "lorentzian", unit = "log2")), sum( log2(1 + abs((P) - (Q)))))
        
})

test_that("distance(method = 'lorentzian') computes the correct distance value using unit = log10.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "lorentzian", unit = "log10")), sum( log10(1 + abs((P) - (Q)))))
        
})


test_that("distance(method = 'intersection') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "intersection")), sum(apply(base::rbind(P,Q),2,min)))
        
})


test_that("distance(method = 'non-intersection') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "non-intersection")), 1 - sum(apply(base::rbind(P,Q),2,min)))
        
})


test_that("distance(method = 'wavehedges') computes the correct distance value.", {
        
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "wavehedges")), sum(abs(P - Q) / apply(base::rbind(P,Q),2,max)))
        
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
        
        expect_equal(as.vector(philentropy::distance(A, B, method = "wavehedges")), wh(A,B))
        
})




test_that("distance(method = 'czekanowski') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "czekanowski")), sum(abs(P - Q)) / sum(P + Q))
        
})


test_that("distance(method = 'motyka') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "motyka")), sum(apply(base::rbind(P,Q),2,max)) / sum(P + Q))
        
})


test_that("distance(method = 'kulczynski_s') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "kulczynski_s")), 1 / (sum(abs((P) - (Q))) / sum(apply(rbind(P, Q),2,min))))
        
})


test_that("distance(method = 'tanimoto') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "tanimoto")), sum(abs((P) - (Q))) / sum(apply(rbind(P, Q),2,max)))
        
})


test_that("distance(method = 'ruzicka') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "ruzicka")), 1 - (sum(abs((P) - (Q))) / sum(apply(rbind(P, Q),2,max))))
        
})


test_that("distance(method = 'inner_product') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "inner_product")), sum ( (P) * (Q) ))
        
})


test_that("distance(method = 'harmonic_mean') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "harmonic_mean")), 2 * sum ( (P) * (Q) / ((P) + (Q)) ))
        
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
        
        expect_equal(as.vector(philentropy::distance(A, B, method = "harmonic_mean")), hm(A,B))        
})


test_that("distance(method = 'cosine') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "cosine")), ((sum ( (P) * (Q) )) / (sqrt(sum((P)^2)) * sqrt(sum((Q)^2)))))
        
})


test_that("distance(method = 'hassebrook') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "hassebrook")), ((sum ( (P) * (Q) )) / (sum((P)^2) + sum((Q)^2) - ((sum ( (P) * (Q) ))))))
        
})



test_that("distance(method = 'jaccard') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "jaccard")), 1 - ((sum ( (P) * (Q) )) / (sum((P)^2) + sum((Q)^2) - ((sum ( (P) * (Q) ))))))
        
})


test_that("distance(method = 'dice') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "dice")), 1 - (2 * (sum ( (P) * (Q) )) / (sum((P)^2) + sum((Q)^2) )))
        
})


test_that("distance(method = 'fidelity') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "fidelity")), sum(sqrt(P * Q)))
        
})


test_that("distance(method = 'bhattacharyya') computes the correct distance value using unit = log.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "bhattacharyya")), -log(sum(sqrt(P * Q))))
        
})

test_that("distance(method = 'bhattacharyya') computes the correct distance value using unit = log2.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "bhattacharyya", unit = "log2")), -log2(sum(sqrt(P * Q))))
        
})

test_that("distance(method = 'bhattacharyya') computes the correct distance value using unit = log10.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "bhattacharyya", unit = "log10")), -log10(sum(sqrt(P * Q))))
        
})


test_that("distance(method = 'hellinger') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "hellinger")), 2 * sqrt(1 - sum(sqrt(P * Q))))
        
})



test_that("distance(method = 'matusita') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "matusita")), sqrt(sum((sqrt(P) - sqrt(Q))^2)))
        
})


test_that("distance(method = 'squared_chord') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "squared_chord")), sum((sqrt(P) - sqrt(Q))^2))
        
})


test_that("distance(method = 'squared_euclidean') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "squared_euclidean")), sum(((P) - (Q))^2))
        
})


test_that("distance(method = 'pearson') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "pearson")), sum(((P) - (Q))^2 / (Q)))
        
})


test_that("distance(method = 'neyman') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "neyman")), sum(((P) - (Q))^2 / (P)))
        
})


test_that("distance(method = 'squared_chi') computes the correct distance value.", {
        
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "squared_chi")), sum(((P) - (Q))^2 / ((P) + (Q))))
        
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
        
        expect_equal(as.vector(philentropy::distance(A, B, method = "squared_chi")), sqchisq(A,B))
        
})



test_that("distance(method = 'prob_symm') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "prob_symm")), 2 * sum(((P) - (Q))^2 / ((P) + (Q))))
        
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
        
        expect_equal(as.vector(philentropy::distance(A, B, method = "prob_symm")), probsymmchisq(A,B))
        
})

test_that("distance(method = 'divergence') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "divergence")), 2 * sum(((P) - (Q))^2 / ((P) + (Q))^2))
        
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
        
        expect_equal(as.vector(philentropy::distance(A, B, method = "divergence")), div(A,B))
})

test_that("distance(method = 'clark') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "clark")), sqrt(sum((abs((P) - (Q)) / ((P) + (Q)))^2)))
        
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
        
        expect_equal(as.vector(philentropy::distance(A, B, method = "clark")), clark(A,B))
        
})


test_that("distance(method = 'additive_symm') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "additive_symm")), sum((((P) - (Q))^2 * ((P) + (Q))) / ((P) * (Q))))
        
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
        
        expect_equal(as.vector(philentropy::distance(A, B, method = "additive_symm")), add(A,B))
})

test_that("distance(method = 'kullback-leibler') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "kullback-leibler")), sum((P) * log((P) / (Q))))
        
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
        
        expect_equal(as.vector(philentropy::distance(A, B, method = "kullback-leibler")), kl(A,B))
        
})




test_that("distance(method = 'jeffreys') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "jeffreys")), sum(((P) - (Q)) * log((P) / (Q))))
        
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
        expect_equal(as.vector(philentropy::distance(A, B, method = "jeffreys")), jeff(A,B))
})



test_that("distance(method = 'k_divergence') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "k_divergence")), sum((P) * log(2 * (P) / ((P) + (Q)))))
        
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
        expect_equal(as.vector(philentropy::distance(A, B, method = "k_divergence")), kdiv(A,B))
})


test_that("distance(method = 'topsoe') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "topsoe")), sum(((P) * log(2 * (P) / ((P) + (Q)))) + ((Q) * log(2 * (Q) / ((P) + (Q))))))
        
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
        expect_equal(as.vector(philentropy::distance(A, B, method = "topsoe")), topsoe(A,B))
        
        
})

test_that("distance(method = 'jensen-shannon') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "jensen-shannon")), 0.5 * ((sum((P) * log((2 * (P)) / ((P) + (Q)))))  +  (sum((Q) * log((2 * (Q)) / ((P) + (Q)))))))
        
})

test_that("distance(method = 'jensen-shannon') computes the correct distance value in case input probability vectors store 0 values at the same position causing 0 * log(0) computation
.", {
        
        A <- c(0,0.25,0.25,0,0.25,0.25)
        B <- c(0,0,0.25,0.25,0.25,0.25)
        
        js <- function(x,y){
                
                dist <- vector(mode = "numeric", length = 1)
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
        expect_equal(as.vector(philentropy::distance(A, B, method = "jensen-shannon")), js(A,B))
})


test_that("distance(method = 'jensen_difference') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "jensen_difference")), sum(((((P) * log((P))) + ((Q) * log((Q)))) / 2 ) - (((P) + (Q)) / 2) * log(((P) + (Q)) / 2)))
        
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
        expect_equal(as.vector(philentropy::distance(A, B, method = "jensen_difference")), js.diff(A,B))
       
})

test_that("distance(method = 'taneja') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "taneja")), sum(((P + Q) / 2) * log((P+Q) / (2 * sqrt(P*Q)))) )
        
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
        expect_equal(as.vector(philentropy::distance(A, B, method = "taneja")), taneja(A,B))
        
})

test_that("distance(method = 'kumar-johnson') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "kumar-johnson")), sum(((P^2 - Q^2)^2 / (2 * (P*Q)^1.5))) )
        
})



test_that("distance(method = 'avg') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(P, Q, method = "avg")), (sum( abs(P-Q) ) + max( abs(P-Q))) / 2 )
        
})

