#' @title Distance between Probability Density Functions
#' @description This functions computes the distance/dissimilarity between two probability density functions.
#' @param x a numeric vector (probability density function).
#' @param y a numeric vector (probability density function).
#' @param method a character string specifying the distance measure that shall be computed.
#' @author Hajk-Georg Drost
#' @details The following distance measures are implemented in this function:
#' @examples
#' 
#' distance(1:10, 20:29, method = "euclidean")
#' 
#' @export

distance <- function(x,y, method = "euclidean", p = NULL){
        
        if(is.na(x) || is.na(y)){
                stop("Your vector stores NA values...")
        }
        
        # result distance
        dist <- NA_real_
        
        if(method == "euclidean"){
                
                dist <- euclidean(x,y)
                
        }
        
        
        if(method == "cityblock"){
                
                dist <- cityblock(x,y)
                
        }
           
        
        if(method == "minkowski"){
                
                if(!is.null(p)){
                        
                        dist <- minkowski(x,y, p)
                        
                } else {
                        
                        stop("Please specify n for the Minkowski distance!")
                }
        }
        
        
        if(method == "chebyshev"){
                
                dist <- chebyshev(x,y)
                
        }
        
        
        if(method == "sorensen"){
                
                dist <- sorensen(x,y)
        }
        
        
        names(dist) <- method
        
        return(dist)
}






