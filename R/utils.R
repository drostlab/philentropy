#' @useDynLib philentropy
#' @importFrom Rcpp sourceCpp
NULL



#' @title Check the validity of the input probability distribution
valid.distr <- function(x){
        
        
        if(any(is.na(x))){
                
                stop("Your input vector includes NA values", call. = FALSE)
        }
        
        else if(!all(dplyr::between(x,0,1))){
                
                stop("Your probability values are not between: [0,1].", call. = FALSE)
        }
        
        else if(abs(1 - sum(x)) > 1E-5) {
                
                stop("Your probability distribution does not sum to 1.", call. = FALSE)
        } else{
                
           return(TRUE)
        
        }
}

#' @title Get method names for \code{distance}
#' @description This function returns the names of the methods that can be applied
#' to compute distances between probability density functions using the \code{\link{distance}} 
#' function.
#' @author Hajk-Georg Drost
#' @examples
#' 
#' getDistMethods()
#' 
#' @export

getDistMethods <- function(){
        
        distance.names <- c("euclidean", "manhattan", "minkowski", "chebyshev",
                            "sorensen", "gower", "soergel", "kulczynski_d",
                            "canberra", "lorentzian", "intersection", "non-intersection",
                            "wavehedges", "czekanowski", "motyka","kulczynski_s",
                            "tanimoto", "ruzicka","inner_product","harmonic_mean",
                            "cosine", "hassebrook", "jaccard", "dice","fidelity","bhattacharyya",
                            "hellinger", "matusita", "squared_chord","squared_euclidean","pearson",
                            "neyman")
        
        return(distance.names)
}
