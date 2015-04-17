#  Part of the philentropy package
#
#  Copyright (C) 2015 Hajk-Georg Drost
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/




#' @title Distances between Probability Density Functions
#' @description This functions computes the distance/dissimilarity between two probability density functions.
#' @param x a numeric \code{data.frame} or \code{matrix} (storing probability vectors) or a count \code{data.frame} or \code{matrix} (if \code{est.prob = TRUE}).
#' @param method a character string indicating whether the distance measure that should be computed.
#' @param p power of the Minkowski distance.
#' @param test.na a boolean value indicating whether input vectors should be tested for \code{NA} values.
#' @param unit a character string specifying the logarithm unit that should be used to compute distances that depend on log computations.
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
#' \item Clark : \eqn{d = sqrt ( \sum (| P_i - Q_i | / (P_i + Q_i))^2 )}
#' \item Additive Symmetric \eqn{\Chi}^2 : \eqn{d = \sum ( ((P_i - Q_i)^2 * (P_i + Q_i)) / (P_i * Q_i) ) }
#' }
#' 
#' \item Shannon's entropy family
#' \itemize{
#' \item Kullback-Leibler : \eqn{d = \sum P_i * log(P_i / Q_i)}
#' \item Jeffreys : \eqn{d = \sum (P_i - Q_i) * log(P_i / Q_i)}
#' \item K divergence : \eqn{d = \sum P_i * log(2 * P_i / P_i + Q_i)}
#' \item Topsoe : \eqn{d = \sum ( P_i * log(2 * P_i / P_i + Q_i) ) + ( Q_i * log(2 * Q_i / P_i + Q_i) )}
#' \item Jensen-Shannon :  \eqn{d = 0.5 * ( \sum P_i * log(2 * P_i / P_i + Q_i) + \sum Q_i * log(2 * Q_i / P_i + Q_i))} 
#' \item Jensen difference : \eqn{d = \sum ( (P_i * log(P_i) + Q_i * log(Q_i) / 2) - (P_i + Q_i / 2) * log(P_i + Q_i / 2) )}
#' }
#' 
#' \item Combinations
#' \itemize{
#' \item Taneja : \eqn{d = \sum ( P_i + Q_i / 2) * log( P_i + Q_i / ( 2 * sqrt( P_i * Q_i)) )}
#' \item Kumar-Johnson : \eqn{d = \sum (P_i^2 - Q_i^2)^2 / 2 * (P_i * Q_i)^1.5}
#' \item Avg(L_1, L_n) : \eqn{d = \sum | P_i - Q_i| + max{ | P_i - Q_i |} / 2}
#' }
#' 
#' }
#' @examples
#' 
#' # receive a list of implemented probability distance measures
#' getDistMethods()
#' 
#' # Simple Examples
#' 
#' ## compute the euclidean distance between two probability vectors
#' distance(rbind(1:10/sum(1:10), 20:29/sum(20:29)), method = "euclidean")
#' 
#' ## compute the euclidean distance between all pairwise comparisons of probability vectors
#' ProbMatrix <- rbind(1:10/sum(1:10), 20:29/sum(20:29),30:39/sum(30:39))
#' distance(ProbMatrix, method = "euclidean")
#' 
#' # compute distance matrix without testing for NA values in the input matrix
#' distance(ProbMatrix, method = "euclidean", test.na = FALSE)
#' @note According to the reference in some distance measure computations invalid computations can
#' occur when dealing with 0 probabilities. 
#' 
#' In these cases the convention is treated as follows:
#' 
#' \itemize{
#' \item division by zero - case \code{0/0}: when the divisor and dividend become zero, \code{0/0} is treated as \code{0}. 
#' \item division by zero - case \code{n/0}: when only the divisor becomes \code{0}, the corresponsning \code{0} is replaced by a small \eqn{\epsilon = 0.00001}.
#' \item log of zero - case \code{0 * log(0)}: is treated as \code{0}.
#' \item log of zero - case \code{log(0)}: zero is replaced by a small \eqn{\epsilon = 0.00001}.
#' }
#' 
#' @return 
#' The following results are returned depending on the dimension of \code{x}:
#' 
#' \itemize{
#' \item in case \code{nrow(x)} = 2 : a single distance value.
#' \item in case \code{nrow(x)} > 2 : a distance \code{matrix} storing distance values for all pairwise probability vector comparisons.  
#' }
#' 
#' @export

distance <- function(x ,
                     method  = "euclidean", 
                     p       = NULL, 
                     test.na = TRUE, 
                     unit    = "log"){
        
        nrows <- NA_integer_
        nrows <- nrow(x)
        
        if(nrows < 2)
                stop("Your input matrix stores only one probability or count vector.")
        
        if(!is.element(method,getDistMethods()))
                stop("Method '",method,"' is not implemented in this function. Please consult getDistMethods().")
        
        if(!is.numeric(x))
                stop("Non numeric values cannot be used to compute distances..")
                
        if(!is.element(unit,c("log","log2","log10")))
                stop("You can only choose units: log, log2, or log10.")
        
        if(!is.element(class(x),c("data.frame","matrix")))
                stop("x should be either a numeric matrix or a numeric data.frame")
        
        if(class(x) == "data.frame")
                x <- as.matrix(x)
        
        
        # although validation would be great, it cost a lot of computation time
        # for large comparisons between multiple distributions
        # here a smarter (faster) way to validate distributions needs to be implemented
#         valid.distr(x)
#         valid.distr(y)
        
        if(nrows == 2){
                # result distance
                dist <- NA_real_
        } else {
                dist <- matrix(NA_real_, nrows, nrows)
        }
        
        
        if(method == "euclidean"){
                
                if(nrows == 2){
                        dist <- euclidean(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,euclidean,test.na)
                }
                
        }
        
        
        if(method == "manhattan"){
                
                if(nrows == 2){
                        dist <- manhattan(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,manhattan,test.na)
                }
                
        }
           
        
        if(method == "minkowski"){
                
                if(!is.null(p)){
                        
                        if(nrows == 2){
                                dist <- minkowski(x[1, ],x[2, ],p,test.na)
                        } else {
                                dist <- DistMatrixMinkowski(x,p,test.na)
                        }
                        
                } else {
                        
                        stop("Please specify p for the Minkowski distance!")
                }
        }
        
        
        if(method == "chebyshev"){
                
                if(nrows == 2){
                        dist <- chebyshev(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,chebyshev,test.na)
                }
        }
        
        
        if(method == "sorensen"){
                
                if(nrows == 2){
                        dist <- sorensen(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,sorensen,test.na)
                }
        }
        
        
        if(method == "gower"){
                
                if(nrows == 2){
                        dist <- gower(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,gower,test.na)
                }
        }
        
        if(method == "soergel"){
                
                if(nrows == 2){
                        dist <- soergel(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,soergel,test.na)
                }
        }
        
        if(method == "kulczynski_d"){
                
                if(nrows == 2){
                        dist <- kulczynski_d(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,kulczynski_d,test.na)
                }
        }
        
        if(method == "canberra"){
                
                if(nrows == 2){
                        dist <- canberra(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,canberra,test.na)
                }
        }

        if(method == "lorentzian"){
                
                if(nrows == 2){
                        dist <- lorentzian(x[1, ],x[2, ],test.na, unit)
                } else {
                        dist <- DistMatrixWithUnit(x,lorentzian,test.na, unit)
                }
        }
        
        if(method == "intersection"){
                
                if(nrows == 2){
                        dist <- intersection_dist(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,intersection_dist,test.na)
                }
        }
        
        if(method == "non-intersection"){
                
                if(nrows == 2){
                        dist <- 1.0 - intersection_dist(x[1, ],x[2, ],test.na)
                } else {
                        dist <- 1.0 - DistMatrixWithoutUnit(x,intersection_dist,test.na)
                }
        }
        
        if(method == "wavehedges"){
                
                if(nrows == 2){
                        dist <- wave_hedges(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,wave_hedges,test.na)
                }
        }
        
        if(method == "czekanowski"){
                
                if(nrows == 2){
                        dist <- czekanowski(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,czekanowski,test.na)
                }
        }
        
        if(method == "motyka"){
                
                if(nrows == 2){
                        dist <- motyka(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,motyka,test.na)
                }
        }
        
        if(method == "kulczynski_s"){
                
                if(nrows == 2){
                        dist <- 1.0 / kulczynski_d(x[1, ],x[2, ],test.na)
                } else {
                        dist <- 1.0 / DistMatrixWithoutUnit(x,kulczynski_d,test.na)
                }
        }
        
        if(method == "tanimoto"){
                
                if(nrows == 2){
                        dist <- tanimoto(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,tanimoto,test.na)
                } 
        }
        
        if(method == "ruzicka"){
                
                if(nrows == 2){
                        dist <- ruzicka(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,ruzicka,test.na)
                }
        }
        
        if(method == "inner_product"){
                
                if(nrows == 2){
                        dist <- inner_product(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,inner_product,test.na)
                }
        }
        
        if(method == "harmonic_mean"){
                
                if(nrows == 2){
                        dist <- harmonic_mean_dist(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,harmonic_mean_dist,test.na)
                }
        }
        
        if(method == "cosine"){
                
                if(nrows == 2){
                        dist <- cosine_dist(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,cosine_dist,test.na)
                }       
        }
        
        if(method == "hassebrook"){
                
                if(nrows == 2){
                        dist <- kumar_hassebrook(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,kumar_hassebrook,test.na)
                } 
        }
        
        if(method == "jaccard"){
                
                if(nrows == 2){
                        dist <- jaccard(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,jaccard,test.na)
                }  
        }
        
        if(method == "dice"){
                
                if(nrows == 2){
                        dist <- dice_dist(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,dice_dist,test.na)
                } 
        }
        
        if(method == "fidelity"){
                
                if(nrows == 2){
                        dist <- fidelity(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,fidelity,test.na)
                } 
        }
        
        if(method == "bhattacharyya"){
                
                if(nrows == 2){
                        dist <- bhattacharyya(x[1, ],x[2, ],test.na,unit)
                } else {
                        dist <- DistMatrixWithUnit(x,bhattacharyya,test.na,unit)
                }
        }
        
        if(method == "hellinger"){
                
                if(nrows == 2){
                        dist <- hellinger(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,hellinger,test.na)
                }     
        }
        
        if(method == "matusita"){
                
                if(nrows == 2){
                        dist <- matusita(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,matusita,test.na)
                }  
        }
        
        if(method == "squared_chord"){
                
                if(nrows == 2){
                        dist <- squared_chord(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,squared_chord,test.na)
                }
        }

        if(method == "squared_euclidean"){
                
                if(nrows == 2){
                        dist <- squared_euclidean(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,squared_euclidean,test.na)
                }       
        }
        
        if(method == "pearson"){
                
                if(nrows == 2){
                        dist <- pearson_chi_sq(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,pearson_chi_sq,test.na)
                }
        }
        
        if(method == "neyman"){
                
                if(nrows == 2){
                        dist <- neyman_chi_sq(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,neyman_chi_sq,test.na)
                } 
        }
        
        if(method == "squared_chi"){
                
                if(nrows == 2){
                        dist <- squared_chi_sq(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,squared_chi_sq,test.na)
                }         
        }
        
        if(method == "prob_symm"){
                
                if(nrows == 2){
                        dist <- prob_symm_chi_sq(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,prob_symm_chi_sq,test.na)
                }      
        }
        
        if(method == "divergence"){
                
                if(nrows == 2){
                        dist <- divergence_sq(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,divergence_sq,test.na)
                }    
        }
        
        if(method == "clark"){
                
                if(nrows == 2){
                        dist <- clark_sq(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,clark_sq,test.na)
                }    
        }
        
        if(method == "additive_symm"){
                
                if(nrows == 2){
                        dist <- additive_symm_chi_sq(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,additive_symm_chi_sq,test.na)
                }
        }
        
        if(method == "kullback-leibler"){
                
                if(nrows == 2){
                        dist <- kullback_leibler_distance(x[1, ],x[2, ],test.na,unit)
                } else {
                        dist <- DistMatrixWithUnit(x,kullback_leibler_distance,test.na,unit)
                }
        }
        
        if(method == "jeffreys"){
                
                if(nrows == 2){
                        dist <- jeffreys(x[1, ],x[2, ],test.na,unit)
                } else {
                        dist <- DistMatrixWithUnit(x,jeffreys,test.na,unit)
                }
        }
        
        if(method == "k_divergence"){
                
                if(nrows == 2){
                        dist <- k_divergence(x[1, ],x[2, ],test.na,unit)
                } else {
                        dist <- DistMatrixWithUnit(x,k_divergence,test.na,unit)
                }
        }
        
        if(method == "topsoe"){
                
                if(nrows == 2){
                        dist <- topsoe(x[1, ],x[2, ],test.na,unit)
                } else {
                        dist <- DistMatrixWithUnit(x,topsoe,test.na,unit)
                }
        }
        
        if(method == "jensen-shannon"){
                
                if(nrows == 2){
                        dist <- jensen_shannon(x[1, ],x[2, ],test.na,unit)
                } else {
                        dist <- DistMatrixWithUnit(x,jensen_shannon,test.na,unit)
                }
        }
        
        if(method == "jensen_difference"){
                
                if(nrows == 2){
                        dist <- jensen_difference(x[1, ],x[2, ],test.na,unit)
                } else {
                        dist <- DistMatrixWithUnit(x,jensen_difference,test.na,unit)
                }
        }
        
        if(method == "taneja"){
                
                if(nrows == 2){
                        dist <- taneja(x[1, ],x[2, ],test.na,unit)
                } else {
                        dist <- DistMatrixWithUnit(x,taneja,test.na,unit)
                }
        }
        
        if(method == "kumar-johnson"){
                
                if(nrows == 2){
                        dist <- kumar_johnson(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,kumar_johnson,test.na)
                }
        }
        
        if(method == "avg"){
                
                if(nrows == 2){
                        dist <- avg(x[1, ],x[2, ],test.na)
                } else {
                        dist <- DistMatrixWithoutUnit(x,avg,test.na)
                }
        }
        
        if(nrows == 2){
                names(dist) <- method
        } else {
                colnames(dist) <- paste0("pvec.",1:nrows)
                rownames(dist) <- paste0("pvec.",1:nrows)
        }
        
        return(dist)
}






