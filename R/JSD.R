#' @title Jensen-Shannon Divergence 
#' @description This function computes the Jensen-Shannon Divergence of two probability distributions P and Q with equal weights.
#' @param P a probability distribution P.
#' @param Q a probability distribution Q.
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
#' JSD.Example <- JSD(P,Q)
#' 
#' # reading a standard PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#' 
#' Prob <- Probability(PhyloExpressionSetExample)
#' 
#' # compute the JSD comparing the PS probability distributions of developmental stages 1 and 5
#' JSD( P = Prob[ , 1],
#'      Q = Prob[ , 5] )
#' @seealso
#' \code{\link{KL}}, \code{\link{KL.Matrix}}, \code{\link{E}}, \code{\link{gJSD}}
#' @export

JSD <- function(P,Q){
        
        valid.distr(P)
        valid.distr(Q)
        
        if(!(length(P) == length(Q)))
                stop("Your input vectors have different lengths..")
        
        return(JensonShannonDivergenceCpp(as.vector(P),as.vector(Q)))
        
}