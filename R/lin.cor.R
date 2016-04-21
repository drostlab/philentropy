#' @title Linear Correlation
#' @description This function computed the linear correlation between two vectors or a correlation matrix for an input matrix.
#' 
#' The following methods to compute linear correlations are implemented in this function:
#'  
#' @param x a numeric \code{vector}, \code{matrix}, or \code{data.frame}.
#' @param y a numeric \code{vector} that should be correlated with \code{x}.
#' @param method the method to compute the linear correlation between \code{x} and \code{y}.
#' @param test.na a boolean value indicating whether input data should be checked for \code{NA} values.
#' @author Hajk-Georg Drost
#' @details
#'  
#'    \itemize{
#'   \item \code{method = "pearson"} : Pearson's correlation coefficient (centred).
#'   \item \code{method = "pearson2"} : Pearson's uncentred correlation coefficient.
#'   \item \code{method = "sq_pearson"} . Squared Pearson's correlation coefficient.
#'   \item \code{method = "kendall"} : Kendall's correlation coefficient.
#'   \item \code{method = "spearman"} : Spearman's correlation coefficient.
#'   } 
#'   
#'  Further Details:
#'   
#'  \itemize{
#'  \item \emph{Pearson's correlation coefficient (centred)} : 
#'  } 
#' @export 
lin.cor <- function(x,y = NULL, method = "pearson", test.na = FALSE){
        
        if(is.null(y)){
                if(!is.element(class(x), c("matrix","data.frame")))
                   stop("x should be either a data.frame or matrix object.")
                   
                cor.coef <- matrix(NA_real_,ncol(x),ncol(x))
        } else {
                if((!is.vector(x)) && (!is.vector(y)))
                        stop("x and y should be vectors.")
                
                cor.coef <- vector("numeric",1)
        }
        
        if(method == "pearson"){
                
                if(class(x) == "matrix"){
                        cor.coef <-  DistMatrixWithoutUnitMAT(x,pearson_corr_centred,test.na)
                }
                else if ((is.vector(x)) && (is.vector(y))){
                        cor.coef <- pearson_corr_centred(x,y,test.na)
                }
        }
        
        if(method == "pearson2"){
                
                if(class(x) == "matrix"){
                        cor.coef <-  DistMatrixWithoutUnitMAT(x,pearson_corr_uncentred,test.na)
                } else if ((is.vector(x)) && (is.vector(y))){
                        cor.coef <-  pearson_corr_uncentred(x,y,test.na)
                }
        }
        
        if(method == "sq_pearson"){
                
                if(class(x) == "matrix"){
                        cor.coef <-  DistMatrixWithoutUnitMAT(x,squared_pearson_corr,test.na)
                } else if ((is.vector(x)) && (is.vector(y))){
                        cor.coef <-  squared_pearson_corr(x,y,test.na)
                }
        }
        
        if(class(cor.coef) == "matrix"){
                colnames(cor.coef) <- paste0("vec.",1:ncol(cor.coef))
                rownames(cor.coef) <- paste0("vec.",1:ncol(cor.coef))       
        } else if (is.vector(cor.coef)){
                names(cor.coef) <- method
        }
        
        return (cor.coef)
}





