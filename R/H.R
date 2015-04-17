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


#' Shannon's Entropy \eqn{H(X)} 
#' 
#' Compute the Shannon's Entropy \eqn{H(X) = - \sum P(X) * log2(P(X))} based on a
#' given probability vector \eqn{P(X)}.
#' 
#' This function might be useful to fastly compute Shannon's Entropy for any
#' given probability vector.
#' 
#' @param x a numeric probability vector \eqn{P(X)} for which
#' Shannon's Entropy \eqn{H(X)} shall be computed.
#' @param unit a character string specifying the logarithm unit that shall be used to compute distances that depend on log computations.
#' @return a numeric value representing Shannon's Entropy in bit.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{JE}}, \code{\link{CE}}, \code{\link{KL}}, \code{\link{JSD}}, \code{\link{gJSD}}
#' @references Shannon, Claude E. 1948. "A Mathematical Theory of
#' Communication". \emph{Bell System Technical Journal} \bold{27} (3): 379-423.
#' @examples
#' 
#' H(1:10/sum(1:10))
#' 
#' @export

H <- function(x, unit = "log2"){
        
        if(anyNA(x))
                stop("x includes NA values... ")
        
        valid.distr(x)
        return(Ecpp(as.vector(x), unit))
        
}
