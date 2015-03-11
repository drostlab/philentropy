#' Generate a Probability Matrix
#' 
#' This function computes the probability of each
#' Phylostratum/Divergence Stratum for each developmental stage \code{s}. The
#' probability of a given phylostratum \code{ps} is defined as the sum of expression
#' levels in stage \code{s} of genes corresponding to \code{ps} divided by the overall sum of
#' expression levels in stage \code{s}.
#' 
#' \deqn{P_s(PS=ps) = \sum ( e_is * \delta(ps_i = ps) ) / \sum e_is}
#' 
#' Where \eqn{e_is} denotes the expression level of gene \eqn{i} in stage
#' \bold{\eqn{s}}, \bold{\eqn{\delta}} denotes the \bold{Kronecker delta}, and
#' \bold{\eqn{ps_i}} denotes the corresponding phylostratum of gene \eqn{i}.
#' 
#' The goal of this function is to model the discrete PS or DS values as probabilities.
#' 
#' @param ExpressionSet a standard PhyloExpressionSet or
#' DivergenceExpressionSet object.
#' @return a numeric matrix storing the probability distributions of PS or DS
#' for each developmental stage s.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{Phylotranscriptomics_Probability}},
#' \code{\link{Phylotranscriptomics_JointProbability}},
#' \code{\link{JointProbability}}, \code{\link{Phylotranscriptomics_Entropy}},
#' \code{\link{Entropy}}
#' @examples
#' 
#' # example data sets
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#' 
#' # generate the Phylostratum Probability Matrix for the given PhyloExpressionSet
#' Probability(PhyloExpressionSetExample)
#' 
#' # generate the Divergence Stratum Probability Matrix for the given DivergenceExpressionSet
#' Probability(DivergenceExpressionSetExample)
#' 
#' @export

Probability <- function(ExpressionSet){
        
        is.ExpressionSet(ExpressionSet)
        
        nCols <- dim(ExpressionSet)[2]
        nRows <- dim(ExpressionSet)[1]
        PhylostratumVector <- vector(mode="numeric",length = nRows)
        PhylostratumVector <- ExpressionSet[ , 1]
        nPS <- as.numeric(max(PhylostratumVector))
        ProbMatrix <- matrix(NA_real_, nrow = nPS, ncol = nCols-2)
                
        ProbMatrix <- Phylotranscriptomics_Probability(as.matrix(ExpressionSet[ , 3:nCols]), PhylostratumVector, nPS)
                
        rownames(ProbMatrix) <- 1:nPS
        colnames(ProbMatrix) <- names(ExpressionSet)[3:nCols]
        return(ProbMatrix)
        
}




