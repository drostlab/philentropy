#' Kullback–Leibler Divergence
#' 
#' This function computes the Kullback–Leibler divergence of two probability
#' distributions P and Q.
#' 
#' \deqn{KL(P||Q) = \sum P(P) * log2(P(P) / P(Q)) = H(P,Q) - H(P)}
#' 
#' where H(P,Q) denotes the \code{\link{JointEntropy}} of the probability
#' distributions P and Q and H(P) denotes the \code{\link{Entropy}} of
#' probability distribution P. In case P = Q then KL(P,Q) = 0 and in case P !=
#' Q then KL(P,Q) > 0.
#' 
#' The KL divergence is a \textbf{non-symmetric measure} of the directed divergence
#' between two probability distributions P and Q. It only fulfills the
#' \emph{positivity} property of a \emph{distance metric}.
#' 
#' Because of the relation KL(P||Q) = H(P,Q) - H(P), the Kullback–Leibler
#' divergence of two probability distributions P and Q is also named
#' \emph{Cross Entropy} of two probability distributions P and Q.
#' @param P a probability distribution P.
#' @param Q a probability distribution Q.
#' @return The Kullback–Leibler divergence of P and Q.
#' @author Hajk-Georg Drost
#' @seealso
#' \code{\link{KL.Matrix}}, \code{\link{Entropy}}, \code{\link{JointEntropy}}
#' @references Cover Thomas M. and Thomas Joy A. 2006. "Elements of Information
#' Theory". \emph{John Wiley & Sons}.
#' @examples
#' 
#' # a general example: comparing a normal distribution with an exponential distribution
#' P <- pnorm(1:10, 5 , 1)
#' Q <- pexp(1:10, 1)
#' KLD <- KL(Q,P)
#' 
#' # a phylotranscriptomics example:
#' # compare the probability distribution of developmental 
#' # stage 1 and 5 of a given PhyloExpressionSet
#' 
#' # read standard phylotranscriptomics data
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#' 
#' pKLD <- KL( P = Probability(PhyloExpressionSetExample)[ , 1],
#'             Q = Probability(PhyloExpressionSetExample)[ , 5] )
#' 
#' 
#'@export
 
KL <- function(P,Q){
        
        valid.distr(P)
        valid.distr(Q)
        
        if(!(length(P) == length(Q)))
                stop("Your input vectors have different lengths..")
        
        return(CrossEntropy(as.vector(P),as.vector(Q)))
        
}