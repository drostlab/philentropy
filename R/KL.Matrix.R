#' Kullback-Leibler Divergence of a Probability Matrix
#' 
#' This function computes the cross entropy (Kullback-Leibler Divergence) of
#' all combinations of stage specific probability distributions.
#' 
#' @param x a probability matrix.
#' @return a matrix of pairwise KL divergences between all combinations of possible comparisons.
#' @author Hajk-Georg Drost
#' @seealso
#' \code{\link{gJSD}}, \code{\link{KL}}, \code{\link{Probability}}
#' @references Cover Thomas M. and Thomas Joy A. 2006. "Elements of Information
#' Theory". \emph{John Wiley & Sons}.
#' @examples
#' 
#' # read standard phylotranscriptomics data
#' data(PhyloExpressionSetExample)
#' 
#' # compute the KL matrix of a given PhyoExpressionSet
#' KLMatrix <- KL.Matrix(Probability(PhyloExpressionSetExample))
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

        if(is.matrix(x) == TRUE){
                
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



