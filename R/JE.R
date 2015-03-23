#' Shannon's Joint-Entropy \eqn{H(X,Y)} 
#' 
#' This funciton computes Shannon's Joint-Entropy \eqn{H(X,Y) = - \sum \sum P(X,Y) *
#' log2(P(X,Y))} based on a given joint-probability vector \eqn{P(X,Y)}.
#' 
#' @param x a numeric joint-probability vector \eqn{P(X,Y)} for
#' which Shannon's Joint-Entropy \eqn{H(X,Y)} shall be computed.
#' @return a numeric value representing Shannon's Joint-Entropy in bit.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{H}}, \code{\link{CE}}, \code{\link{KL}}, \code{\link{KL.Matrix}}, \code{\link{JSD}}, \code{\link{gJSD}}, \code{\link{distance}}
#' @references Shannon, Claude E. 1948. "A Mathematical Theory of
#' Communication". \emph{Bell System Technical Journal} \bold{27} (3): 379-423.
#' @examples
#' 
#' JE(1:100/sum(1:100))
#' 
#' @export

JE <- function(x){
        
       valid.distr(x)
       return(JEcpp(as.vector(x)))
       
}
