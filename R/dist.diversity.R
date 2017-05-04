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




#' Distance Diversity between Probability Density Functions
#' 
#' This function computes all distance values between two probability density functions that are available in \code{\link{getDistMethods}}
#' and returns a vector storing the corresponding distance measures. This vector is \emph{named distance diversity vector}.
#' 
#' @param x a numeric \code{data.frame} or \code{matrix} (storing probability vectors) or a numeric \code{data.frame} or \code{matrix} storing counts (if \code{est.prob} is specified).
#' @param p power of the Minkowski distance.   
#' @author Hajk-Georg Drost
#' @examples
#' 
#' dist.diversity(rbind(1:10/sum(1:10), 20:29/sum(20:29)))
#' 
#' @export
        
dist.diversity <- function(x, p = 2){
        
        distMethods <- vector(mode = "character")
        nMethods <- NA_integer_
        distMethods <- getDistMethods()
        nMethods <- length(distMethods)
        distDiversityVec <- vector(mode = "numeric", length = nMethods)
        
        distDiversityVec <- unlist(lapply(distMethods, function(method) distance(x, method = method, p = as.double(p))))
        names(distDiversityVec) <- distMethods
        
        return(distDiversityVec)

}