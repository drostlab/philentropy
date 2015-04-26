#' @title Linear Correlation
#' @description This function computed the linear correlation between two vectors or a correlation matrix for an input matrix.
#' 
#' The following methods to compute linear correlations are implemented in this function:
#'  
#'  @param x a numeric \code{vector}, \code{matrix}, or \code{data.frame}.
#'  @param y a numeric \code{vector} that should be correlated with \code{x}.
#'  @author Hajk-Georg Drost
#'  @details
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
#'    @export 
lin.cor <- function(x,y = NULL, method = "pearson"){
        
        cor.coef <- vector("numeric",1)
                
        if(method == "pearson"){
                cor.coef <-  pearson_corr_centred(x,y)
        }
        
        if(method == "pearson2"){
                cor.coef <-  pearson_corr_noncentred(x,y)
        }
        
        if(method == "sq_pearson"){
                cor.coef <-  squared_pearson_corr(x,y)
        }
        
        return (cor.coef)
}





