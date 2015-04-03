#' Distance Diversity between Probability Density Functions
#' 
#' This function computes all distance values between two probability density functions that are available in \code{\link{getDistMethods}}
#' and returns a vector storing the corresponding distance measures. This vector is \emph{named distance diversity vector}.
#' 
#' @param x a numeric vector (probability density function).
#' @param y a numeric vector (probability density function).
#' @param p power of the Minkowski distance.   
#' 
#' @author Hajk-Georg Drost
#' @examples
#' 
#' dist.diversity(1:10/sum(1:10), 20:29/sum(20:29))
#' 
#' @export
        
dist.diversity <- function(x,y, p = 2){
        
        distMethods <- vector(mode = "character")
        nMethods <- NA_integer_
        distMethods <- getDistMethods()
        nMethods <- length(distMethods)
        distDiversityVec <- vector(mode = "numeric", length = nMethods)
        
        distDiversityVec <- sapply(distMethods, function(method) distance(x = x, y = y, method = method, p = p))
        names(distDiversityVec) <- distMethods
        
        return(distDiversityVec)

}