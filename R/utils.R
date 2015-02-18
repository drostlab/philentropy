#' @useDynLib phylentropy
#' @importFrom Rcpp sourceCpp
NULL

#' @title Testing the validaty of a PhyloExpressionSet or DivergenceExpressionSet standard
#' @description This function tests whether a given ExpressionSet follows the pre-defined PhyloExpressionSet or DivergenceExpressionSet standard.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet that shall be tested for format validity.
#' @author Hajk-Georg Drost
#' @examples
#' 
#' # read example PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#' 
#' is.ExpressionSet(PhyloExpressionSetExample)
#' 
#' @export
is.ExpressionSet <- function(ExpressionSet){
        
        ncols <- dim(ExpressionSet)[2]
        
        d.f_bool <- is.data.frame(ExpressionSet)
        age.vector_bool <- is.numeric(ExpressionSet[ , 1])
        gene.vector_bool <- ifelse(is.factor(ExpressionSet[ , 2]),is.character(levels(ExpressionSet[ , 2])),is.character(ExpressionSet[ , 2]))
        expression.matrix_bool <- all(sapply(ExpressionSet[ , 3:ncols], is.numeric))
        any.NA.values_bool <- !any(is.na(ExpressionSet))
        
        
        if(all(c(d.f_bool,age.vector_bool,gene.vector_bool,expression.matrix_bool,any.NA.values_bool))){
                return(TRUE)
        }
        
        else{
                stop("The present input object does not fulfill the ExpressionSet standard.")
        }
}

