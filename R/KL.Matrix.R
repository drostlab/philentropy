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


#' Kullback-Leibler Divergence of a Probability Matrix
#' 
#' This function computes the cross entropy (Kullback-Leibler Divergence) of
#' all combinations of stage specific probability distributions.
#' 
#' @param x a probability matrix.
#' @return a matrix of pairwise KL divergences between all combinations of possible comparisons.
#' @author Hajk-Georg Drost
#' @seealso
#' \code{\link{gJSD}}, \code{\link{KL}}, \code{\link{H}}, \code{\link{CE}}, \code{\link{JE}}
#' @references Cover Thomas M. and Thomas Joy A. 2006. "Elements of Information
#' Theory". \emph{John Wiley & Sons}.
#' @examples
#' 
#' Prob <- cbind(1:10/sum(1:10), 20:29/sum(20:29), 30:39/sum(30:39))
#' 
#' # compute the KL matrix of a given probability matrix
#' KLMatrix <- KL.Matrix(Prob)
#' 
#' # plot a heatmap of the corresponding KL matrix
#' heatmap(KLMatrix)
#' 
#' @export

KL.Matrix <- function(x){
        
        # check for ditribution validity
        apply(x,2,valid.distr)
        
        nCols <- dim(x)[2]
        KLMatrix <- matrix(NA_real_,nCols,nCols)

        if(is.matrix(x)){
                
                for(i in 1:nCols){
                        for(j in 1:nCols){
                                KLMatrix[i,j] <- CrossEntropy(as.vector(x[ , i]),as.vector(x[ , j]))
                        } 
                }
                
                rownames(KLMatrix) <- colnames(x)
                colnames(KLMatrix) <- colnames(x)
                
                return(KLMatrix)
                
        } else{
                stop("Please use a numeric matrix as input.")
        }
}



