#' @title Jensen-Shannon Divergence 
#' @description This function computes the Jensen-Shannon Divergence of two probability distributions P and Q with equal weights.
#' @param P a probability distribution P.
#' @param Q a probability distribution Q.
#' @param test.na a boolean value specifying whether input vectors shall be tested for NA values.
#' @param unit a character string specifying the logarithm unit that shall be used to compute distances that depend on log computations.
#' @return The Jensen-Shannon divergence  of P and Q.
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
#' \deqn{JSP(P || Q) = 0.5 * (KL(P || R) + KL(Q || R))}
#' 
#' where \eqn{R = 0.5 * (P + Q)} denotes the mid-point of the probability
#' vectors P and Q, and KL(P || R), KL(Q || R) denote the Kullback-Leibler
#' Divergence of P and R, as well as Q and R.
#' 
#' General properties of the Jensen-Shannon Divergence:
#' 
#' 1) JSD is non-negative.
#' 
#' 2) JSD is a symmetric measure JSD(P || Q) = JSD(Q || P).
#' 
#' 3) JSD = 0, if and only if P = Q.
#' 
#' @references Lin J. 1991. "Divergence Measures Based on the Shannon Entropy".
#' IEEE Transactions on Information Theory. (33) 1: 145-151.
#' 
#' Endres M. and Schindelin J. E. 2003. "A new metric for probability
#' distributions". IEEE Trans. on Info. Thy. (49) 3: 1858-1860.
#' @examples
#' 
#' P <- 1:10/sum(1:10)
#' Q <- 20:29/sum(20:29)
#' 
#' JSD.Example <- JSD(P,Q)
#' 
#' @seealso
#' \code{\link{KL}}, \code{\link{H}}, \code{\link{CE}}, \code{\link{gJSD}}, \code{\link{distance}}
#' @export

JSD <- function(P,Q, test.na = TRUE, unit = "log2"){
                
        return( distance(x      = P,
                         y      = Q, 
                        method  = "jensen-shannon", 
                        test.na = test.na, 
                        unit    = unit) )
        
}





