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

#' Shannon's Mutual Information \eqn{I(X,Y)}
#' 
#' Compute Shannon's Mutual Information based on the identity \eqn{I(X,Y) =
#' H(X) + H(Y) - H(X,Y)} based on a given joint-probability vector \eqn{P(X,Y)}
#' and probability vectors \eqn{P(X)} and \eqn{P(Y)}.
#' 
#' This function might be useful to fastly compute Shannon's Mutual Information
#' for any given joint-probability vector and probability vectors.
#' 
#' @param x a numeric probability vector \eqn{P(X)}.
#' @param y a numeric probability vector \eqn{P(Y)}.
#' @param xy a numeric joint-probability vector \eqn{P(X,Y)}.
#' @param unit a character string specifying the logarithm unit that shall be used to compute distances that depend on log computations.
#' @return Shannon's Mutual Information in bit.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{H}}, \code{\link{JE}}, \code{\link{CE}}
#' @references Shannon, Claude E. 1948. "A Mathematical Theory of
#' Communication". \emph{Bell System Technical Journal} \bold{27} (3): 379-423.
#' @examples
#' 
#' MI( x = 1:10/sum(1:10), y = 20:29/sum(20:29), xy = 1:10/sum(1:10) )
#' 
#' @export

MI <- function(x,y,xy, unit = "log2"){
        
        valid.distr(x)
        valid.distr(y)
        valid.distr(xy)
        
        return(MIcpp(as.vector(x),as.vector(y),as.vector(xy), unit))
        
}






