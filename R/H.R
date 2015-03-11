#' Shannon's Entropy \eqn{H(X)} 
#' 
#' Compute the Shannon's Entropy \eqn{H(X) = - \sum P(X) * log2(P(X))} based on a
#' given probability vector \eqn{P(X)}.
#' 
#' This function might be useful to fastly compute Shannon's Entropy for any
#' given probability vector.
#' 
#' @param x a numeric probability vector \eqn{P(X)} for which
#' Shannon's Entropy \eqn{H(X)} shall be computed.
#' @return a numeric value representing Shannon's Entropy in bit.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{CE}}, \code{\link{KL}}, \code{\link{JSD}}, \code{\link{gJSD}}, \code{\link{Entropy}}
#' @references Shannon, Claude E. 1948. "A Mathematical Theory of
#' Communication". \emph{Bell System Technical Journal} \bold{27} (3): 379-423.
#' @examples
#' 
#' # read a standard PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#' 
#' # compute the Entropy based on the probability P(PS) of developmental stage 1
#' EntropyValue <- H(Probability(PhyloExpressionSetExample)[ , 1])
#' 
#' # or
#' 
#' # compute the phylotranscriptomics Entropy profile analogous to the 'Entopy' function
#' EntropyProfile <- apply(Probability(PhyloExpressionSetExample),2,H)
#' 
#' # and compare it with the 'Entropy' function
#' EntropyProfile2 <- Entropy(PhyloExpressionSetExample)
#' 
#' @export

H <- function(x){
        
        valid.distr(x)
        return(Ecpp(as.vector(x)))
        
}
