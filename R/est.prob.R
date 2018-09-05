#' @title Estimate Probability Vectors From Count Vectors
#' @description This function takes a numeric count vector and returns estimated
#' probabilities of the corresponding counts.
#' 
#' The following probability estimation methods are implemented in this function:
#' 
#' \itemize{
#' \item \code{method = "empirical"} : generates the relative frequency of the data \code{x/sum(x)}.
#' \item
#' \item
#' }
#' 
#' @param x a numeric vector storing count values.
#' @param method a character string specifying the estimation method tht should be used to estimate probabilities from input counts.
#' @author Hajk-Georg Drost
#' @return a numeric probability vector.
#' @examples
#' 
#' # generate a count vector
#' x <- runif(100)
#' 
#' # generate a probability vector from corresponding counts
#' # method = "empirical"
#' x.prob <- estimate.probability(x, method = "empirical")
#' 
#' @export


estimate.probability <- function(x, method = "empirical"){
        
        if(!is.element(method,c("empirical")))
                stop("Please choose a valid probability estimation method.")
        
        if(method == "empirical"){
                # fastest implementation for relative frequency 
                return(x/sum(x))
        }
        
}