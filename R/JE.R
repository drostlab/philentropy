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


#' Shannon's Joint-Entropy \eqn{H(X,Y)} 
#' 
#' This funciton computes Shannon's Joint-Entropy \eqn{H(X,Y) = - \sum \sum P(X,Y) *
#' log2(P(X,Y))} based on a given joint-probability vector \eqn{P(X,Y)}.
#' 
#' @param x a numeric joint-probability vector \eqn{P(X,Y)} for
#' which Shannon's Joint-Entropy \eqn{H(X,Y)} shall be computed.
#' @param unit a character string specifying the logarithm unit that shall be used to compute distances that depend on log computations.
#' @return a numeric value representing Shannon's Joint-Entropy in bit.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{H}}, \code{\link{CE}}, \code{\link{KL}}, \code{\link{JSD}}, \code{\link{gJSD}}, \code{\link{distance}}
#' @references Shannon, Claude E. 1948. "A Mathematical Theory of
#' Communication". \emph{Bell System Technical Journal} \bold{27} (3): 379-423.
#' @examples
#' 
#' JE(1:100/sum(1:100))
#' 
#' @export

JE <- function(x, unit = "log2"){
        
       valid.distr(x)
       return(JEcpp(as.vector(x), unit))
       
}
