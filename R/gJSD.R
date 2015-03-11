#' @title Generalized Jensen-Shannon Divergence 
#' @description This function computes the Generalized Jensen-Shannon Divergence of a probability matrix with equal weights.
#' @param P a probability distribution P.
#' @param Q a probability distribution Q.
#' @return The Jensen-Shannon divergence  of P and Q.
#' @author Hajk-Georg Drost
#' @details 
#' 
#' Function to compute the Generalized Jensen-Shannon Divergence
#' 
#' @examples
#' 
#' Prob <- Probability(PhyloExpressionSetExample)
#' 
#' # compute the Generalized JSD comparing the PS probability matrix
#' gJSD(Prob)
#' 
#' @seealso
#' \code{\link{KL}}, \code{\link{KL.Matrix}}, \code{\link{E}}, \code{\link{JSD}}
#' @export

gJSD <- function(ProbabilityMatrix){
        
        # check for ditribution validity
        apply(ProbabilityMatrix,2,valid.distr)
        
        nDistributions <- dim(ProbabilityMatrix)[2]
        nElements <- dim(ProbabilityMatrix)[1]
        # defining the weights for the generalized Jensen-Shannon Divergence
        weights <- NA_real_
        g.JSD <- NA_real_ 
        weights <- rep(1/nDistributions, nDistributions)
        weightedProbabilityMatrix <- matrix(NA_real_, nrow = nElements, ncol = nDistributions)
        
        
        if(length(weights) == nDistributions){
                
                for(i in 1:nDistributions){
                        
                        weightedProbabilityMatrix[ ,i] <- weights[i] * ProbabilityMatrix[ ,i]
                        
                }
                
                g.JSD <- H(rowSums(weightedProbabilityMatrix)) - sum((weights * apply(ProbabilityMatrix,2,H)))
                
                return(g.JSD)
                
        } else{
                stop("The weight vector and probability distribution have different lengths!")
        }
}







