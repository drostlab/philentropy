#' Shannon's Joint-Entropy \eqn{H(X,Y)} 
#' 
#' This funciton computes Shannon's Joint-Entropy \eqn{H(X,Y) = - \sum \sum P(X,Y) *
#' log2(P(X,Y))} based on a given joint-probability vector \eqn{P(X,Y)}.
#' 
#' @param x a numeric joint-probability vector \eqn{P(X,Y)} for
#' which Shannon's Joint-Entropy \eqn{H(X,Y)} shall be computed.
#' @return a numeric value representing Shannon's Joint-Entropy in bit.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{H}}, \code{\link{KL}}, \code{\link{JSD}}, \code{\link{distance}}
#' @references Shannon, Claude E. 1948. "A Mathematical Theory of
#' Communication". \emph{Bell System Technical Journal} \bold{27} (3): 379-423.
#' @examples
#' 
#' # reading a standard PhyloExpressionSet and standard DivergenceExpressionSet
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#' 
#' # compute the Joint-Entropy based on the joint-probability P(PS,DS) of developmental stage 1
#' JE(JointProbability(PhyloExpressionSetExample,DivergenceExpressionSetExample)[ , 1])
#' 
#' # or
#' 
#' # compute the phylotranscriptomics Entropy profile analogous to the 'JointEntopy' function
#' apply(JointProbability(PhyloExpressionSetExample,DivergenceExpressionSetExample),2,JE)
#' 
#' # and compare it with the 'Entropy' function:
#' 
#' # generating a PhyloDivergenceExpressionSet object
#' Ex <- getPhyloDivergenceExpressionSet(PhyloExpressionSetExample,DivergenceExpressionSetExample)
#' 
#' JointEntropy(Ex)
#' 
#' @export

JE <- function(x){
        
       valid.distr(x)
       return(JEcpp(as.vector(x)))
       
}
