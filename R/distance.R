#' @title Distance between Probability Density Functions
#' @description This functions computes the distance/dissimilarity between two probability density functions.
#' @param x a numeric vector (probability density function).
#' @param y a numeric vector (probability density function).
#' @param method a character string specifying the distance measure that shall be computed.
#' @author Hajk-Georg Drost
#' @details The following distance measures are implemented in this function:
#' 
#' \itemize{
#' \item L_p Minkowski family
#' \itemize{
#' \item Euclidean : \eqn{d = sqrt( \sum | P_i - Q_i |^2)}
#' \item City block (manhatten) :
#' \item Minkowski :
#' \item Chebyshev : 
#' }
#' 
#' \item L_1 family
#' \itemize{
#' \item Sorensen :
#' \item Gower : 
#' }
#' 
#' \item Shannon's entropy family
#' \itemize{
#' \item Kullback-Leibler : \eqn{KL(P || Q) = \sum P(P) * log2(P(P) / P(Q)) = H(P,Q) - H(P)}
#' \item Jensen-Shannon : \eqn{JSP(P || Q) = 0.5 * (KL(P || R) + KL(Q || R))}, where \eqn{R = 0.5 * (P + Q)} denotes the mid-point of the probability
#' vectors P and Q
#' }
#' }
#' @examples
#' 
#' distance(1:10, 20:29, method = "euclidean")
#' 
#' @export

distance <- function(x,y, method = "euclidean", p = NULL){
        
        if(is.na(x) || is.na(y)){
                stop("Your vector stores NA values...")
        }
        
        if(!all(is.numeric(x),is.numeric(y))){
                stop("Non numeric values cannot be used to compute distances..")
                
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
        
        
        
        if(method == "kullback-leibler"){
                
                dist <- KL(x,y)
        }
        
        
        if(method == "jensen-shannon"){
                
                dist <- JSD(x,y)
        }
        
        names(dist) <- method
        
        return(dist)
}






