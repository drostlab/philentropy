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



#' @useDynLib philentropy, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function(libpath) {
  library.dynam.unload("philentropy", libpath)
}

# @title Check the validity of input probability distributions
valid.distr <- function(x){
        
        sapply(x, function(vec) if(!dplyr::between(vec,0,1)) stop("Your probability values are not between: [0,1].", .call = FALSE))
        
        
        if (anyNA(x)){
                stop("Your input vector includes NA values", call. = FALSE)
        }
        
        else if (sum_rcpp(x) > 1.0000001) {
                
                stop("Your probability distribution does not sum to 1.", call. = FALSE)
        } 
}

#' @title Get method names for \code{distance}
#' @description This function returns the names of the methods that can be applied
#' to compute distances between probability density functions using the \code{\link{distance}} 
#' function.
#' @author Hajk-Georg Drost
#' @examples
#' 
#' getDistMethods()
#' 
#' @export

getDistMethods <- function(){
        
        distance.names <- vector(mode = "character", length = 46)
        distance.names <- c("euclidean", "manhattan", "minkowski", "chebyshev",
                            "sorensen", "gower", "soergel", "kulczynski_d",
                            "canberra", "lorentzian", "intersection", "non-intersection",
                            "wavehedges", "czekanowski", "motyka","kulczynski_s",
                            "tanimoto", "ruzicka","inner_product","harmonic_mean",
                            "cosine", "hassebrook", "jaccard", "dice","fidelity","bhattacharyya",
                            "hellinger", "matusita", "squared_chord","squared_euclidean","pearson",
                            "neyman", "squared_chi", "prob_symm", "divergence","clark",
                            "additive_symm","kullback-leibler","jeffreys","k_divergence",
                            "topsoe","jensen-shannon", "jensen_difference","taneja",
                            "kumar-johnson","avg")
        
        return(distance.names)
}
