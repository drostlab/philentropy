#' @title Generalized Jensen-Shannon Divergence 
#' @description This function computes the Generalized Jensen-Shannon Divergence of a probability matrix with equal weights.
#' @param x a probability matrix.
#' @return The Jensen-Shannon divergence between all possible combinations of comparisons.
#' @author Hajk-Georg Drost
#' @details 
#' 
#' Function to compute the Generalized Jensen-Shannon Divergence
#' 
#' @examples
#' 
#' Prob <- cbind(1:10/sum(1:10), 20:29/sum(20:29), 30:39/sum(30:39))
#' 
#' # compute the Generalized JSD comparing the PS probability matrix
#' gJSD(Prob)
#' 
#' @seealso
#' \code{\link{KL}}, \code{\link{KL.Matrix}}, \code{\link{H}}, \code{\link{JSD}},
#'  \code{\link{CE}}, \code{\link{JE}}
#' @export

gJSD <- function(x){
        
        # check for ditribution validity
        apply(x,2,valid.distr)
        
        nDistributions <- dim(x)[2]
        nElements <- dim(x)[1]
        # defining the weights for the generalized Jensen-Shannon Divergence
        weights <- NA_real_
        g.JSD <- NA_real_ 
        weights <- rep(1/nDistributions, nDistributions)
        weightedProbabilityMatrix <- matrix(NA_real_, nrow = nElements, ncol = nDistributions)
        
        
        if(length(weights) == nDistributions){
                
                for(i in 1:nDistributions){
                        
                        weightedProbabilityMatrix[ ,i] <- weights[i] * x[ ,i]
                        
                }
                
                g.JSD <- H(rowSums(weightedProbabilityMatrix)) - sum((weights * apply(x,2,H)))
                
                return(g.JSD)
                
        } else{
                stop("The weight vector and probability distribution have different lengths!")
        }
}







