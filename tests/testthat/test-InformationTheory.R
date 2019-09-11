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


context("Test implementation of gJSD() ...")

test_that("Numeric computation of gJSD() : ", {

        skip_on_cran()
        
        P <- 1:10/sum(1:10)
        Q <- 20:29/sum(20:29)
        R <- 30:39/sum(30:39)
        
        x <- cbind(P,Q,R)
        w <- rep(1 / ncol(x),ncol(x))
        wpm <- w * x
        gjsd <- H(rowSums(wpm)) - sum((w * apply(x,2,H)))
        
        expect_equal(gJSD(rbind(P,Q,R)), gjsd)
        
})


test_that("gJSD() internally changes a data.frame to a matrix", {
        
        skip_on_cran()
        
        P <- 1:10/sum(1:10)
        Q <- 20:29/sum(20:29)
        R <- 30:39/sum(30:39)
        
        x <- cbind(P,Q,R)
        w <- rep(1 / ncol(x),ncol(x))
        wpm <- w * x
        gjsd <- H(rowSums(wpm)) - sum((w * apply(x,2,H)))
        
        expect_equal(gJSD(data.frame(rbind(P,Q,R))), gjsd)
        
})

test_that("gJSD() checks for transposed matrix column sums > 1.", {
        
expect_equal(gJSD(matrix(c(1, 1, 0, 0), nrow = 2)), 0)
expect_equal(gJSD(matrix(c(1, 0, 0, 1), nrow = 2)), 1)
        
})
