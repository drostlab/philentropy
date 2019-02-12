#' @title Jensen-Shannon Divergence 
#' @description This function computes a distance matrix or distance value based on the Jensen-Shannon Divergence with equal weights.
#' @param x a numeric \code{data.frame} or \code{matrix} (storing probability vectors) or a numeric \code{data.frame} or \code{matrix} storing counts (if \code{est.prob = TRUE}). See \code{\link{distance}} for details.
#' @param test.na a boolean value specifying whether input vectors shall be tested for NA values.
#' @param unit a character string specifying the logarithm unit that shall be used to compute distances that depend on log computations.
#' @param est.prob method to estimate probabilities from input count vectors such as non-probability vectors. Default: \code{est.prob = NULL}. Options are:
#' \itemize{
#' \item \code{est.prob = "empirical"}: The relative frequencies of each vector are computed internally. For example an input matrix \code{rbind(1:10, 11:20)} will be transformed to a probability vector \code{rbind(1:10 / sum(1:10), 11:20 / sum(11:20))}
#' }
#' @return a distance value or matrix based on JSD computations.
#' @author Hajk-Georg Drost
#' @details 
#' 
#' Function to compute the Jensen-Shannon Divergence JSD(P || Q) between two
#' probability distributions P and Q with equal weights \eqn{\pi_1} =
#' \eqn{\pi_2} = \eqn{1/2}.
#' 
#' The Jensen-Shannon Divergence JSD(P || Q) between two probability
#' distributions P and Q is defined as:
#' 
#' \deqn{JSD(P || Q) = 0.5 * (KL(P || R) + KL(Q || R))}
#' 
#' where \eqn{R = 0.5 * (P + Q)} denotes the mid-point of the probability
#' vectors P and Q, and KL(P || R), KL(Q || R) denote the Kullback-Leibler
#' Divergence of P and R, as well as Q and R.
#' 
#' \strong{General properties of the Jensen-Shannon Divergence:}
#' 
#' \itemize{
#' \item \code{1)} JSD is non-negative.
#' \item \code{2)} JSD is a symmetric measure JSD(P || Q) = JSD(Q || P).
#' \item \code{3)} JSD = 0, if and only if P = Q.
#' }
#' @references Lin J. 1991. "Divergence Measures Based on the Shannon Entropy".
#' IEEE Transactions on Information Theory. (33) 1: 145-151.
#' 
#' Endres M. and Schindelin J. E. 2003. "A new metric for probability
#' distributions". IEEE Trans. on Info. Thy. (49) 3: 1858-1860.
#' @examples
#' # Jensen-Shannon Divergence between P and Q
#' P <- 1:10/sum(1:10)
#' Q <- 20:29/sum(20:29)
#' x <- rbind(P,Q)
#' JSD(x)
#' 
#' # Jensen-Shannon Divergence between P and Q using different log bases
#' JSD(x, unit = "log2") # Default
#' JSD(x, unit = "log")
#' JSD(x, unit = "log10")
#' 
#' # Jensen-Shannon Divergence Divergence between count vectors P.count and Q.count
#' P.count <- 1:10
#' Q.count <- 20:29
#' x.count <- rbind(P.count,Q.count)
#' JSD(x.count, est.prob = "empirical")
#' 
#' # Example: Distance Matrix using JSD-Distance
#' 
#' Prob <- rbind(1:10/sum(1:10), 20:29/sum(20:29), 30:39/sum(30:39))
#'
#' # compute the KL matrix of a given probability matrix
#' JSDMatrix <- JSD(Prob)
#'
#' # plot a heatmap of the corresponding JSD matrix
#' heatmap(JSDMatrix)
#' 
#' @seealso
#' \code{\link{KL}}, \code{\link{H}}, \code{\link{CE}}, \code{\link{gJSD}}, \code{\link{distance}}
#' @export

JSD <- function(x, test.na = TRUE, unit = "log2", est.prob = NULL){
  if (!is.matrix(x))
    stop("Please provide a matrix as input, e.g. with x <- rbind(vector1, vector2).", call. = FALSE)
  
        return( distance(x       = x, 
                        method   = "jensen-shannon", 
                        test.na  = test.na, 
                        unit     = unit,
                        est.prob = est.prob) )
        
}





