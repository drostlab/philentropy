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
#' Default: \code{weights} = \code{NULL}; in this case all distributions are weighted equally (= uniform distribution of weights).
#' In case users wish to specify non-uniform weights for e.g. 3 distributions, they
#' can specify the argument \code{weights = c(0.5, 0.25, 0.25)}. This notation
#' denotes that \code{vec1} is weighted by \code{0.5}, \code{vec2} is weighted by \code{0.25}, and \code{vec3} is weighted by \code{0.25} as well.
#' @param est.prob method to estimate probabilities from input count vectors such as non-probability vectors. Default: \code{est.prob = NULL}. Options are:
#' \itemize{
#' \item \code{est.prob = "empirical"}: The relative frequencies of each vector are computed internally. For example an input matrix \code{rbind(1:10, 11:20)} will be transformed to a probability vector \code{rbind(1:10 / sum(1:10), 11:20 / sum(11:20))}
#' }
#' @return The Jensen-Shannon divergence between all possible combinations of comparisons.
#' @author Hajk-Georg Drost
#' @details
#'
#' Function to compute the Generalized Jensen-Shannon Divergence
#'
#'\eqn{JSD_{\pi_1,...,\pi_n}(P_1, ..., P_n) = H(\sum_{i = 1}^n \pi_i * P_i) - \sum_{i = 1}^n \pi_i*H(P_i)}
#'
#'where \eqn{\pi_1,...,\pi_n} denote the weights selected for the probability vectors \code{P_1,...,P_n} and \code{H(P_i)} denotes the Shannon Entropy of probability vector \code{P_i}.
#' @examples
#' # define input probability matrix
#' Prob <- rbind(1:10/sum(1:10), 20:29/sum(20:29), 30:39/sum(30:39))
#'
#' # compute the Generalized JSD comparing the PS probability matrix
#' gJSD(Prob)
#' 
#' # Generalized Jensen-Shannon Divergence between three vectors using different log bases
#' gJSD(Prob, unit = "log2") # Default
#' gJSD(Prob, unit = "log")
#' gJSD(Prob, unit = "log10")
#' 
#' # Jensen-Shannon Divergence Divergence between count vectors P.count and Q.count
#' P.count <- 1:10
#' Q.count <- 20:29
#' R.count <- 30:39
#' x.count <- rbind(P.count, Q.count, R.count)
#' gJSD(x.count, est.prob = "empirical")
#'
#' @seealso
#' \code{\link{KL}}, \code{\link{H}}, \code{\link{JSD}},
#'  \code{\link{CE}}, \code{\link{JE}}
#' @export

gJSD <- function(x, unit = "log2", weights = NULL, est.prob = NULL) {
        if (is.null(weights))
                message(
                        "No weights were specified ('weights = NULL'), thus equal weights for all distributions will be calculated and applied."
                )

        # transpose input data
        x <- t(as.matrix(x))
        
        message(
                "Metric: 'gJSD'; unit = '",
                unit,
                "'; comparing: ",
                ncol(x),
                " vectors (v1, ... , v",
                ncol(x),
                ")."
        )
        
        if (!is.null(est.prob)) {
                x <- apply(x, 2, estimate.probability, method = est.prob) 
        }
        # check for distribution validity
        apply(x, 2, valid.distr)
        
        nDistributions <- ncol(x)
        nElements <- nrow(x)
        
        if (is.null(weights)) {
                # defining the weights for the generalized Jensen-Shannon Divergence
                # -> weights are equally distributed
                weights <-
                        vector(mode = "numeric", length = nDistributions)
                g.JSD <- NA_real_
                weights <- rep(1 / nDistributions, nDistributions)
                
        } else {
                # check for the validity of input weights
                valid.distr(weights)
                
                if (length(weights) != nDistributions)
                        stop(
                                "The length of your input 'weights' vector differs from the number of input distributions..",
                                call. = FALSE
                        )
                
        }
        
        WeightedProbabilityMatrix <-
                matrix(NA_real_, nrow = nElements, ncol = nDistributions)
        
        for (i in seq_len(length(weights)))
                WeightedProbabilityMatrix[, i] <-
                weights[i] * x[, i]
        
        g.JSD <-
                H(rowSums(WeightedProbabilityMatrix), unit = unit) - sum((weights * apply(x, 2, H, unit = unit)))
        
        message("Weights: ",
                paste0("v", 1:length(weights), " = ", weights, collapse = ", "))
        return(g.JSD)
        
}
