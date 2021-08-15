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

#' Kullback-Leibler Divergence
#' 
#' This function computes the Kullback-Leibler divergence of two probability
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
#' Because of the relation KL(P||Q) = H(P,Q) - H(P), the Kullback-Leibler
#' divergence of two probability distributions P and Q is also named
#' \emph{Cross Entropy} of two probability distributions P and Q.
#' @param x a numeric \code{data.frame} or \code{matrix} (storing probability vectors) or a numeric \code{data.frame} or \code{matrix} storing counts (if \code{est.prob = TRUE}). See \code{\link{distance}} for details.
#' @param test.na a boolean value indicating whether input vectors should be tested for NA values.
#' @param unit a character string specifying the logarithm unit that shall be used to compute distances that depend on log computations.
#' @param est.prob method to estimate probabilities from a count vector. Default: est.prob = NULL.
#' @param epsilon a small value to address cases in the KL computation where division by zero occurs. In
#' these cases, x / 0 or 0 / 0 will be replaced by \code{epsilon}. The default is \code{epsilon = 0.00001}.
#' However, we recommend to choose a custom \code{epsilon} value depending on the size of the input vectors,
#' the expected similarity between compared probability density functions and 
#' whether or not many 0 values are present within the compared vectors.
#' As a rough rule of thumb we suggest that when dealing with very large 
#' input vectors which are very similar and contain many \code{0} values,
#' the \code{epsilon} value should be set even smaller (e.g. \code{epsilon = 0.000000001}),
#' whereas when vector sizes are small or distributions very divergent then
#' higher \code{epsilon} values may also be appropriate (e.g. \code{epsilon = 0.01}).
#' Addressing this \code{epsilon} issue is important to avoid cases where distance metrics
#' return negative values which are not defined and only occur due to the
#' technical issues of computing x / 0 or 0 / 0 cases. 
#' @return The Kullback-Leibler divergence of probability vectors.
#' @author Hajk-Georg Drost
#' @seealso
#' \code{\link{H}}, \code{\link{CE}}, \code{\link{JSD}}, \code{\link{gJSD}}, \code{\link{distance}}
#' @references Cover Thomas M. and Thomas Joy A. 2006. Elements of Information
#' Theory. \emph{John Wiley & Sons}.
#' @examples
#'
#' # Kulback-Leibler Divergence between P and Q
#' P <- 1:10/sum(1:10)
#' Q <- 20:29/sum(20:29)
#' x <- rbind(P,Q)
#' KL(x)
#' 
#' # Kulback-Leibler Divergence between P and Q using different log bases
#' KL(x, unit = "log2") # Default
#' KL(x, unit = "log")
#' KL(x, unit = "log10")
#' 
#' # Kulback-Leibler Divergence between count vectors P.count and Q.count
#' P.count <- 1:10
#' Q.count <- 20:29
#' x.count <- rbind(P.count,Q.count)
#' KL(x, est.prob = "empirical")
#' 
#' # Example: Distance Matrix using KL-Distance
#' 
#' Prob <- rbind(1:10/sum(1:10), 20:29/sum(20:29), 30:39/sum(30:39))
#'
#' # compute the KL matrix of a given probability matrix
#' KLMatrix <- KL(Prob)
#'
#' # plot a heatmap of the corresponding KL matrix
#' heatmap(KLMatrix)
#' 
#' 
#'@export
 
KL <- function(x, test.na = TRUE, unit = "log2", est.prob = NULL, epsilon = 0.00001){
        if (!is.matrix(x))
          stop("Please provide a matrix as input, e.g. with x <- rbind(vector1, vector2).", call. = FALSE)
        
        return( distance( x           = x,
                          method      = "kullback-leibler",
                          test.na     = test.na,
                          unit        = unit,
                          est.prob    = est.prob,
                          epsilon     = epsilon) )
        
}
