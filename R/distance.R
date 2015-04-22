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
#' @param x a numeric \code{data.frame} or \code{matrix} (storing probability vectors) or a numeric \code{data.frame} or \code{matrix} storing counts (if \code{est.prob = TRUE}).
#' @param method a character string indicating whether the distance measure that should be computed.
#' @param p power of the Minkowski distance.
#' @param test.na a boolean value indicating whether input vectors should be tested for \code{NA} values. Faster computations if \code{test.na = FALSE}.
#' @param unit a character string specifying the logarithm unit that should be used to compute distances that depend on log computations.
#' @param est.prob method to estimate probabilities from a count vector. Default: \code{est.prob = NULL}.
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
#' # Simple Examples
#' 
#' # receive a list of implemented probability distance measures
#' getDistMethods()
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
                     method   = "euclidean", 
                     p        = NULL, 
                     test.na  = TRUE, 
                     unit     = "log",
                     est.prob = NULL){
        
        nrows <- NA_integer_
        nrows <- nrow(x)
        
        if(nrows < 2)
                stop("Your input matrix stores only one probability or count vector.")
        
        if(!is.element(method,getDistMethods()))
                stop("Method '",method,"' is not implemented in this function. Please consult getDistMethods().")
        
        
        if(!is.element(class(x),c("data.frame","matrix")))
                stop("x should be a data.frame or matrix.")
        
        
        if(!is.null(est.prob)){
                
                if(class(x) == "data.frame"){
                        # estimate probability row-wise
                        x <- t(apply(x,1,estimate.probability, method = est.prob))
                }
                
                if(class(x) == "matrix"){
                        x <- t(apply(x,1,estimate.probability, method = est.prob))  
                }
        }
        
#         if(!is.numeric(x))
#                 stop("Non numeric values cannot be used to compute distances..")
#                 
        if(!is.element(unit,c("log","log2","log10")))
                stop("You can only choose units: log, log2, or log10.")
        
        
        
        # although validation would be great, it cost a lot of computation time
        # for large comparisons between multiple distributions
        # here a smarter (faster) way to validate distributions needs to be implemented
        # check for distribution validity
       # apply(x,1,valid.distr, test.na = test.na)
        
        
        dist <- matrix(NA_real_, nrows, nrows)
        
        if(method == "euclidean"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,euclidean,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,euclidean,test.na)
                
        }
        
        
        else if(method == "manhattan"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,manhattan,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,manhattan,test.na)
                
        }
           
        
        else if(method == "minkowski"){
                
                if(!is.null(p)){
                        
                        if(class(x) == "data.frame")
                                dist <- DistMatrixMinkowskiDF(x,p,test.na)
                        if(class(x) == "matrix")
                                dist <- DistMatrixMinkowskiMAT(x,p,test.na)
                } else {
                        
                        stop("Please specify p for the Minkowski distance!")
                }
        }
        
        
        else if(method == "chebyshev"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,chebyshev,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,chebyshev,test.na)
        }
        
        
        else if(method == "sorensen"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,sorensen,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,sorensen,test.na)
                
        }
        
        
        else if(method == "gower"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,gower,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,gower,test.na)
        }
        
        else if(method == "soergel"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,soergel,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,soergel,test.na)
        }
        
        else if(method == "kulczynski_d"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,kulczynski_d,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,kulczynski_d,test.na)
        }
        
        else if(method == "canberra"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,canberra,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,canberra,test.na)
        }

        else if(method == "lorentzian"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithUnitDF(x,lorentzian,test.na, unit)
                if(class(x) == "matrix")
                        dist <- DistMatrixWithUnitMAT(x,lorentzian,test.na, unit)
        }
        
        else if(method == "intersection"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,intersection_dist,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,intersection_dist,test.na)
                
        }
        
        else if(method == "non-intersection"){
                
                if(class(x) == "data.frame")
                        dist <- 1.0 - DistMatrixWithoutUnitDF(x,intersection_dist,test.na)
                
                if(class(x) == "matrix")
                        dist <- 1.0 - DistMatrixWithoutUnitMAT(x,intersection_dist,test.na)
        }
        
        else if(method == "wavehedges"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,wave_hedges,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,wave_hedges,test.na)
                
        }
        
        else if(method == "czekanowski"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,czekanowski,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,czekanowski,test.na)
        }
        
        else if(method == "motyka"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,motyka,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,motyka,test.na)
        }
        
        else if(method == "kulczynski_s"){
                
                if(class(x) == "data.frame")
                        dist <- 1.0 / DistMatrixWithoutUnitDF(x,kulczynski_d,test.na)
                
                if(class(x) == "matrix")
                        dist <- 1.0 / DistMatrixWithoutUnitMAT(x,kulczynski_d,test.na)
        }
        
        else if(method == "tanimoto"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,tanimoto,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,tanimoto,test.na)
        }
        
        else if(method == "ruzicka"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,ruzicka,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,ruzicka,test.na)
        }
        
        else if(method == "inner_product"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,inner_product,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,inner_product,test.na)
        }
        
        else if(method == "harmonic_mean"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,harmonic_mean_dist,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,harmonic_mean_dist,test.na)
                
        }
        
        else if(method == "cosine"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,cosine_dist,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,cosine_dist,test.na)
                  
        }
        
        else if(method == "hassebrook"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,kumar_hassebrook,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,kumar_hassebrook,test.na)
                
        }
        
        else if(method == "jaccard"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,jaccard,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,jaccard,test.na)

        }
        
        else if(method == "dice"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,dice_dist,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,dice_dist,test.na)
                
        }
        
        else if(method == "fidelity"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,fidelity,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,fidelity,test.na)
        }
        
        else if(method == "bhattacharyya"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithUnitDF(x,bhattacharyya,test.na,unit)
                if(class(x) == "matrix")
                        dist <- DistMatrixWithUnitMAT(x,bhattacharyya,test.na,unit)
        }
        
        else if(method == "hellinger"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,hellinger,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,hellinger,test.na)  
        }
        
        else if(method == "matusita"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,matusita,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,matusita,test.na)
        }
        
        else if(method == "squared_chord"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,squared_chord,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,squared_chord,test.na)
        }

        else if(method == "squared_euclidean"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,squared_euclidean,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,squared_euclidean,test.na)    
        }
        
        else if(method == "pearson"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,pearson_chi_sq,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,pearson_chi_sq,test.na)    
                
        }
        
        else if(method == "neyman"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,neyman_chi_sq,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,neyman_chi_sq,test.na)
                
        }
        
        else if(method == "squared_chi"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,squared_chi_sq,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,squared_chi_sq,test.na)
                        
        }
        
        else if(method == "prob_symm"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,prob_symm_chi_sq,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,prob_symm_chi_sq,test.na)
                
        }
        
        else if(method == "divergence"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,divergence_sq,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,divergence_sq,test.na)
                  
        }
        
        else if(method == "clark"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,clark_sq,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,clark_sq,test.na)
                
        }
        
        else if(method == "additive_symm"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,additive_symm_chi_sq,test.na)
                
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,additive_symm_chi_sq,test.na)
                
        }
        
        else if(method == "kullback-leibler"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithUnitDF(x,kullback_leibler_distance,test.na,unit)
                if(class(x) == "matrix")
                        dist <- DistMatrixWithUnitMAT(x,kullback_leibler_distance,test.na,unit)
        }
        
        else if(method == "jeffreys"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithUnitDF(x,jeffreys,test.na,unit)
                if(class(x) == "matrix")
                        dist <- DistMatrixWithUnitMAT(x,jeffreys,test.na,unit)
                
        }
        
        else if(method == "k_divergence"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithUnitDF(x,k_divergence,test.na,unit)
                if(class(x) == "matrix")
                        dist <- DistMatrixWithUnitMAT(x,k_divergence,test.na,unit)
        }
        
        else if(method == "topsoe"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithUnitDF(x,topsoe,test.na,unit)
                if(class(x) == "matrix")
                        dist <- DistMatrixWithUnitMAT(x,topsoe,test.na,unit)
                
        }
        
        else if(method == "jensen-shannon"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithUnitDF(x,jensen_shannon,test.na,unit)
                if(class(x) == "matrix")
                        dist <- DistMatrixWithUnitMAT(x,jensen_shannon,test.na,unit)
        }
        
        else if(method == "jensen_difference"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithUnitDF(x,jensen_difference,test.na,unit)
                if(class(x) == "matrix")
                        dist <- DistMatrixWithUnitMAT(x,jensen_difference,test.na,unit)
        }
        
        else if(method == "taneja"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithUnitDF(x,taneja,test.na,unit)
                if(class(x) == "matrix")
                        dist <- DistMatrixWithUnitMAT(x,taneja,test.na,unit)
        }
        
        else if(method == "kumar-johnson"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,kumar_johnson,test.na)
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,kumar_johnson,test.na)
        }
        
        else if(method == "avg"){
                
                if(class(x) == "data.frame")
                        dist <- DistMatrixWithoutUnitDF(x,avg,test.na)
                if(class(x) == "matrix")
                        dist <- DistMatrixWithoutUnitMAT(x,avg,test.na)
        }
        
        if(nrows == 2){
                dist <- as.vector(dist[2,1])
                names(dist) <- method
        } else {
                colnames(dist) <- paste0("pvec.",1:nrows)
                rownames(dist) <- paste0("pvec.",1:nrows)
        }
        
        return(dist)
}






