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


context("Test implementation of information theory measures...")

test_that("Numeric computation of gJSD() : ", {
        
        P <- 1:10/sum(1:10)
        Q <- 20:29/sum(20:29)
        R <- 30:39/sum(30:39)
        
        x <- cbind(P,Q,R)
        w <- rep(1 / ncol(x),ncol(x))
        wpm <- w * x
        gjsd <- H(rowSums(wpm)) - sum((w * apply(x,2,H)))
        
        expect_equal(gJSD(cbind(P,Q,R)), gjsd)
        
})


test_that("gJSD() throughs error when distributions of different lengths (a list) is passed as argument", {
        
        
        A <- 1:10/sum(1:10)
        B <- 20:28/sum(20:28)
        C <- 30:39/sum(30:39)
        P <- 1:10/sum(1:10)
        Q <- 20:29/sum(20:29)
        R <- 30:39/sum(30:39)
        
        x <- cbind(P,Q,R)
        w <- rep(1 / ncol(x),ncol(x))
        wpm <- w * x
        gjsd <- H(rowSums(wpm)) - sum((w * apply(x,2,H)))
        
        expect_error(gJSD(list(A,B,C)),"Please enter a numeric probability matrix.")
})


test_that("gJSD() internally changes a data.frame to a matrix", {
        
        P <- 1:10/sum(1:10)
        Q <- 20:29/sum(20:29)
        R <- 30:39/sum(30:39)
        
        x <- cbind(P,Q,R)
        w <- rep(1 / ncol(x),ncol(x))
        wpm <- w * x
        gjsd <- H(rowSums(wpm)) - sum((w * apply(x,2,H)))
        
        expect_equal(gJSD(data.frame(P,Q,R)), gjsd)
        
})





