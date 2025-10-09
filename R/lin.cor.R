#' @title Linear Correlation
#' @description This function computed the linear correlation between two vectors or a correlation matrix for an input matrix.
#' 
#' The following methods to compute linear correlations are implemented in this function:
#'  
#' @param x a numeric \code{vector}, \code{matrix}, or \code{data.frame}.
#' @param y a numeric \code{vector} that should be correlated with \code{x}.
#' @param method the method to compute the linear correlation between \code{x} and \code{y}.
#' @param test.na a boolean value indicating whether input data should be checked for \code{NA} values.
#' @param epsilon a small value to address cases where division by zero occurs. Default is `0.00001`.
#' This is only used for `method = "pearson"`.
#' @param num.threads an integer specifying the number of threads to be used for parallel computations.
#' Default is `NULL`, which uses the value from the `RCPP_PARALLEL_NUM_THREADS` environment variable or `2` if not set.
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
lin.cor <- function(x, y = NULL, method = "pearson", test.na = FALSE, epsilon = 0.00001, num.threads = NULL){
        
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
                
                if(inherits(x, "matrix")){
                        cor.coef <- distance_cpp(x, "pearson", p = NULL, test_na = test.na, unit = "log", epsilon = epsilon, num_threads = num.threads)
                }
                else if ((is.vector(x)) && (is.vector(y))){
                        cor.coef <- dist_one_one(x, y, method = "pearson", testNA = test.na)
                }
        }
        
        if(method == "pearson2"){
                
                if(inherits(x, "matrix")){
                        stop("The 'pearson2' method is not supported for matrix input in the new parallel backend.")
                } else if ((is.vector(x)) && (is.vector(y))){
                        cor.coef <- pearson_corr_uncentred(x, y, test.na)
                }
        }
        
        if(method == "sq_pearson"){
                
                if(inherits(x, "matrix")){
                        cor.coef <- distance_cpp(x, "squared_pearson", p = NULL, test_na = test.na, unit = "log", epsilon = epsilon, num_threads = num.threads)
                } else if ((is.vector(x)) && (is.vector(y))){
                        cor.coef <- dist_one_one(x, y, method = "squared_pearson", testNA = test.na)
                }
        }
        
        if(inherits(cor.coef, "matrix")){
                colnames(cor.coef) <- paste0("vec.",1:ncol(cor.coef))
                rownames(cor.coef) <- paste0("vec.",1:ncol(cor.coef))       
        } else if (is.vector(cor.coef)){
                names(cor.coef) <- method
        }
        
        return (cor.coef)
}
