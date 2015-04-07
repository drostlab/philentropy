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

#' @title Generalized Jensen-Shannon Divergence 
#' @description This function computes the Generalized Jensen-Shannon Divergence of a probability matrix with equal weights.
#' @param x a probability matrix.
#' @return The Jensen-Shannon divergence between all possible combinations of comparisons.
#' @author Hajk-Georg Drost
#' @details 
#' 
#' Function to compute the Generalized Jensen-Shannon Divergence
#' 
#' @examples
#' 
#' Prob <- cbind(1:10/sum(1:10), 20:29/sum(20:29), 30:39/sum(30:39))
#' 
#' # compute the Generalized JSD comparing the PS probability matrix
#' gJSD(Prob)
#' 
#' @seealso
#' \code{\link{KL}}, \code{\link{KL.Matrix}}, \code{\link{H}}, \code{\link{JSD}},
#'  \code{\link{CE}}, \code{\link{JE}}
#' @export

gJSD <- function(x){
        
        if(class(x) == "data.frame")
                x <- as.matrix(x)
        
        if(!(class(x) == "matrix"))
                stop("Please enter a numeric probability matrix.")
        
        # check for ditribution validity
        apply(x,2,valid.distr)
        
        nDistributions <- ncol(x)
        nElements <- nrow(x)
        # defining the weights for the generalized Jensen-Shannon Divergence
        weights <- vector(mode = "numeric", length = nDistributions)
        g.JSD <- NA_real_ 
        weights <- rep(1/nDistributions, nDistributions)
        weightedProbabilityMatrix <- matrix(NA_real_, nrow = nElements, ncol = nDistributions)
        
        weightedProbabilityMatrix <- weights * x
                
        g.JSD <- H(rowSums(weightedProbabilityMatrix)) - sum((weights * apply(x,2,H)))
                
        return(g.JSD)
        
}







