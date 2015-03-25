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
#' \item Manhattan : \eqn{d = \sum | P_i - Q_i |}
#' \item Minkowski : \eqn{d = ( \sum | P_i - Q_i |^p)^1/p}
#' \item Chebyshev : \eqn{d = max | P_i - Q_i |}
#' }
#' 
#' \item L_1 family
#' \itemize{
#' \item Sorensen : \eqn{d = \sum | P_i - Q_i | / \sum (P_i + Q_i)}
#' \item Gower : \eqn{d = 1/d * \sum | P_i - Q_i |}
#' \item Soergel : \eqn{d = \sum | P_i - Q_i | / \sum max(P_i , Q_i)}
#' \item Kulczynski d : \eqn{d = \sum | P_i - Q_i | / \sum min(P_i , Q_i)}
#' \item Canberra : \eqn{d = \sum | P_i - Q_i | / (P_i + Q_i)}
#' \item Lorentzian : \eqn{d = \sum ln(1 + | P_i - Q_i |)}
#' }
#' 
#' 
#' \item Intersection family
#' \itemize{
#' \item Intersection : \eqn{s = \sum min(P_i , Q_i)}
#' \item Non-Intersection : \eqn{d = 1 - \sum min(P_i , Q_i)}
#' \item Wave Hedges : \eqn{d = \sum | P_i - Q_i | / max(P_i , Q_i)}
#' \item Czekanowski : \eqn{d = \sum | P_i - Q_i | / \sum | P_i + Q_i |}
#' \item Motyka : \eqn{d = \sum min(P_i , Q_i) / (P_i + Q_i)}
#' \item Kulczynski s : \eqn{d = 1 / \sum | P_i - Q_i | / \sum min(P_i , Q_i)}
#' \item Tanimoto : \eqn{d = \sum (max(P_i , Q_i) - min(P_i , Q_i)) / \sum max(P_i , Q_i)} ; equivalent to Soergel
#' \item Ruzicka : \eqn{s = \sum min(P_i , Q_i) / \sum max(P_i , Q_i)} ; equivalent to 1 - Tanimoto = 1 - Soergel 
#' }
#' 
#' \item Inner Product family
#' \itemize{
#' \item Inner Product : \eqn{s = \sum P_i * Q_i}
#' \item Harmonic mean : \eqn{s = 2 * \sum (P_i * Q_i) / (P_i + Q_i)}
#' \item Cosine : \eqn{s = \sum (P_i * Q_i) / sqrt(\sum P_i^2) * sqrt(\sum Q_i^2)}
#' \item Kumar-Hassebrook (PCE) : \eqn{s = \sum (P_i * Q_i) / (\sum P_i^2 + \sum Q_i^2 - \sum (P_i * Q_i))}
#' \item Jaccard : \eqn{d = 1 - \sum (P_i * Q_i) / (\sum P_i^2 + \sum Q_i^2 - \sum (P_i * Q_i))} ; equivalent to 1 - Kumar-Hassebrook
#' \item Dice : \eqn{d = \sum (P_i - Q_i)^2 / (\sum P_i^2 + \sum Q_i^2)}
#' }
#' 
#' \item Squared-chord family
#' \itemize{
#' \item Fidelity : \eqn{s = \sum sqrt(P_i * Q_i)}
#' \item Bhattacharyya : \eqn{d = - ln \sum sqrt(P_i * Q_i)}
#' \item Hellinger : \eqn{d = 2 * sqrt( 1 - \sum sqrt(P_i * Q_i))}
#' \item Matusita : \eqn{d = sqrt( 2 - 2 * \sum sqrt(P_i * Q_i))}
#' \item Squared-chord : \eqn{d = \sum ( sqrt(P_i) - sqrt(Q_i) )^2}
#' }
#' 
#' 
#' \item Squared L_2 family (\eqn{\Chi}^2 squared family)
#' \itemize{
#' \item Squared Euclidean : \eqn{d = \sum ( P_i - Q_i )^2}
#' \item Pearson \eqn{\Chi}^2 : \eqn{d = \sum ( (P_i - Q_i )^2 / Q_i )}
#' \item Neyman \eqn{\Chi}^2 : \eqn{d = \sum ( (P_i - Q_i )^2 / P_i )}
#' \item Squared \eqn{\Chi}^2 : \eqn{d = \sum ( (P_i - Q_i )^2 / (P_i + Q_i) )}
#' \item Probabilistic Symmetric \eqn{\Chi}^2 : \eqn{d = 2 *  \sum ( (P_i - Q_i )^2 / (P_i + Q_i) )}
#' \item Divergence : \eqn{\Chi}^2 : \eqn{d = 2 *  \sum ( (P_i - Q_i )^2 / (P_i + Q_i)^2 )}
#' \item Clark :
#' \item Additive Symmetric \eqn{\Chi}^2 : 
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
#' distance(1:10/sum(1:10), 20:29/sum(20:29), method = "euclidean")
#' 
#' @export

distance <- function(x,y, method = "euclidean", p = NULL){
        
        if(!is.element(method,getDistMethods()))
                stop("Method '",method,"' is not implemented in this function. Please consult getDistMethods().")
        
        if(is.na(x) || is.na(y)){
                stop("Your vector stores NA values...")
        }
        
        if(!all(is.numeric(x),is.numeric(y))){
                stop("Non numeric values cannot be used to compute distances..")
                
        }
        
        # although validation would be great, it cost a lot of computation time
        # for large comparisons between multiple distributions
        # here a smarter (faster) way to validate distributions needs to be implemented
        valid.distr(x)
        valid.distr(y)
        
        # result distance
        dist <- NA_real_
        
        if(method == "euclidean"){
                
                dist <- euclidean(x,y)
                
        }
        
        
        if(method == "manhattan"){
                
                dist <- manhattan(x,y)
                
        }
           
        
        if(method == "minkowski"){
                
                if(!is.null(p)){
                        
                        dist <- minkowski(x,y, p)
                        
                } else {
                        
                        stop("Please specify p for the Minkowski distance!")
                }
        }
        
        
        if(method == "chebyshev"){
                
                dist <- chebyshev(x,y)
                
        }
        
        
        if(method == "sorensen"){
                
                dist <- sorensen(x,y)
        }
        
        
        if(method == "gower"){
                
                dist <- gower(x,y)
                
        }
        
        if(method == "soergel"){
                
                dist <- soergel(x,y)
                
        }
        
        if(method == "kulczynski_d"){
                
                dist <- kulczynski_d(x,y)
                
        }
        
        if(method == "canberra"){
                
                dist <- canberra(x,y)
                
        }

        if(method == "lorentzian"){
                
                dist <- lorentzian(x,y)
                
        }
        
        if(method == "intersection"){
                
                dist <- intersection_dist(x,y)
                
        }
        
        if(method == "non-intersection"){
                
                dist <- 1.0 - intersection_dist(x,y)
                
        }
        
        if(method == "wavehedges"){
                
                dist <- wave_hedges(x,y)
                
        }
        
        if(method == "czekanowski"){
                
                dist <- czekanowski(x,y)
                
        }
        
        if(method == "motyka"){
                
                dist <- motyka(x,y)
                
        }
        
        if(method == "kulczynski_s"){
                
                dist <- 1.0 / kulczynski_d(x,y)
                
        }
        
        if(method == "tanimoto"){
                
                dist <- tanimoto(x,y)
                
        }
        
        if(method == "ruzicka"){
                
                dist <- ruzicka(x,y)
                
        }
        
        if(method == "inner_product"){
                
                dist <- inner_product(x,y)
                
        }
        
        if(method == "harmonic_mean"){
                
                dist <- harmonic_mean_dist(x,y)
                
        }
        
        if(method == "cosine"){
                
                dist <- cosine_dist(x,y)
                
        }
        
        if(method == "hassebrook"){
                
                dist <- kumar_hassebrook(x,y)
                
        }
        
        if(method == "jaccard"){
                
                dist <- jaccard(x,y)
                
        }
        
        if(method == "dice"){
                
                dist <- dice_dist(x,y)
                
        }
        
        if(method == "fidelity"){
                
                dist <- fidelity(x,y)
                
        }
        
        if(method == "bhattacharyya"){
                
                dist <- bhattacharyya(x,y)
                
        }
        
        if(method == "hellinger"){
                
                dist <- hellinger(x,y)
                
        }
        
        if(method == "matusita"){
                
                dist <- matusita(x,y)
                
        }
        
        if(method == "squared_chord"){
                
                dist <- squared_chord(x,y)
                
        }

        if(method == "squared_euclidean"){
                
                dist <- squared_euclidean(x,y)
                
        }
        
        if(method == "pearson"){
                
                dist <- pearson_chi_sq(x,y)
                
        }
        
        if(method == "neyman"){
                
                dist <- neyman_chi_sq(x,y)
                
        }
        
        if(method == "squared_chi"){
                
                dist <- squared_chi_sq(x,y)
                
        }
        
        if(method == "prob_symm"){
                
                dist <- prob_symm_chi_sq(x,y)
                
        }
        
        if(method == "divergence"){
                
                dist <- divergence_sq(x,y)
                
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






