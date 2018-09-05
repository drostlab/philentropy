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


#' Shannon's Conditional-Entropy \eqn{H(X | Y)}
#' 
#' Compute Shannon's Conditional-Entropy based on the chain rule \eqn{H(X | Y)
#' = H(X,Y) - H(Y)} based on a given joint-probability vector \eqn{P(X,Y)} and
#' probability vector \eqn{P(Y)}.
#' 
#' This function might be useful to fastly compute Shannon's
#' Conditional-Entropy for any given joint-probability vector and probability
#' vector.
#' 
#' @note Note that the probability vector P(Y) must be the probability
#' distribution of random variable Y ( P(Y) for which H(Y) is computed ) and
#' furthermore used for the chain rule computation of \eqn{H(X | Y) = H(X,Y) -
#' H(Y)}.
#' @param xy a numeric joint-probability vector \eqn{P(X,Y)}
#' for which Shannon's Joint-Entropy \eqn{H(X,Y)} shall be computed.
#' @param y a numeric probability vector \eqn{P(Y)} for which
#' Shannon's Entropy \eqn{H(Y)} (as part of the chain rule) shall be computed.
#' It is important to note that this probability vector must be the probability
#' distribution of random variable Y ( P(Y) for which H(Y) is computed).
#' @param unit a character string specifying the logarithm unit that shall be used to compute distances that depend on log computations.
#' @return Shannon's Conditional-Entropy in bit.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{H}}, \code{\link{JE}}
#' @references Shannon, Claude E. 1948. "A Mathematical Theory of
#' Communication". \emph{Bell System Technical Journal} \bold{27} (3): 379-423.
#' @examples
#'  
#'  CE(1:10/sum(1:10),1:10/sum(1:10))
#' 
#' @export

CE <- function(xy,y, unit = "log2"){
        valid.distr(xy)
        valid.distr(y)
        return(CEcpp(as.vector(xy),as.vector(y), unit))
        
}

