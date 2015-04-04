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

#' Kullback–Leibler Divergence
#' 
#' This function computes the Kullback–Leibler divergence of two probability
#' distributions P and Q.
#' 
#' \deqn{KL(P||Q) = \sum P(P) * log2(P(P) / P(Q)) = H(P,Q) - H(P)}
#' 
#' where H(P,Q) denotes the joint entropy of the probability
#' distributions P and Q and H(P) denotes the entropy of
#' probability distribution P. In case P = Q then KL(P,Q) = 0 and in case P !=
#' Q then KL(P,Q) > 0.
#' 
#' The KL divergence is a non-symmetric measure of the directed divergence
#' between two probability distributions P and Q. It only fulfills the
#' \emph{positivity} property of a \emph{distance metric}.
#' 
#' Because of the relation KL(P||Q) = H(P,Q) - H(P), the Kullback–Leibler
#' divergence of two probability distributions P and Q is also named
#' \emph{Cross Entropy} of two probability distributions P and Q.
#' @param P a probability distribution P.
#' @param Q a probability distribution Q.
#' @return The Kullback–Leibler divergence of P and Q.
#' @author Hajk-Georg Drost
#' @seealso
#' \code{\link{H}}, \code{\link{CE}}, \code{\link{JSD}}
#' @references Cover Thomas M. and Thomas Joy A. 2006. "Elements of Information
#' Theory". \emph{John Wiley & Sons}.
#' @examples
#' 
#' # a general example: comparing a normal distribution with an exponential distribution
#' P <- 1:10/sum(1:10)
#' Q <- 20:29/sum(20:29)
#' KLD <- KL(P,Q)
#' 
#'@export
 
KL <- function(P,Q){
        
        valid.distr(P)
        valid.distr(Q)
        
        if(!(length(P) == length(Q)))
                stop("Your input vectors have different lengths..")
        
        return(CrossEntropy(as.vector(P),as.vector(Q)))
        
}
