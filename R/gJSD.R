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
#' @description This function computes the Generalized Jensen-Shannon Divergence of a probability matrix.
#' @param x a probability matrix.
#' @param unit a character string specifying the logarithm unit that shall be used to compute distances that depend on log computations.
#' @param weights a numeric vector specifying the weights for each distribution in \code{x}. 
#' 
#' Default: \code{weights} = \code{NULL}; in this case all distributions are weighted equally. 
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
#' \code{\link{KL}}, \code{\link{H}}, \code{\link{JSD}},
#'  \code{\link{CE}}, \code{\link{JE}}
#' @export

gJSD <- function(x, unit = "log2", weights = NULL){
        
        if(class(x) == "data.frame")
                x <- as.matrix(x)
        
        if(!(class(x) == "matrix"))
                stop("Please enter a numeric probability matrix.")
        
        # check for distribution validity
        apply(x,2,valid.distr)
        
        nDistributions <- ncol(x)
        nElements <- nrow(x)
        
        if(is.null(weights)){
                
                # defining the weights for the generalized Jensen-Shannon Divergence
                # -> weights are equally distributed
                weights <- vector(mode = "numeric", length = nDistributions)
                g.JSD <- NA_real_ 
                weights <- rep(1/nDistributions, nDistributions)
                
        } else {
                # check for the validity of input weights
                valid.distr(weights)
                
                if(length(weights) != nDistributions)
                        stop("The length of your input 'weights' vector differs from the number of input distributions..")
                
        }
        
        WeightedProbabilityMatrix <- matrix(NA_real_, nrow = nElements, ncol = nDistributions)
        
        WeightedProbabilityMatrix <- weights * x
                
        g.JSD <- H(rowSums(WeightedProbabilityMatrix), unit = unit) - sum((weights * apply(x,2,H, unit = unit)))
                
        return(g.JSD)
        
}







