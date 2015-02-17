#' @title Perform Sampling for Entropy
#' 
#' @description
#' This function receives a standard PhyloExpressionSet or
#' DivergenceExpressionSet, samples the Phylostratum or Divergence Stratum column of the corresponding
#' PhyloExpressionSet or DivergenceExpressionSet and furthermore, computes the
#' Entropy profile of this sampled PhyloExpressionSet or
#' DivergenceExpressionSet.
#' 
#' The resulting sample matrix consists of \code{N x S} elements, where S denotes the
#' total number of developmental stages stored in the corresponding
#' PhyloExpressionSet or DivergenceExpressionSet and \code{N} the number of permutations. This \code{N x S} sample matrix can
#' be used to perform further test statistics (sampling test) to quantify the
#' statistical significance of any observed Entropy pattern.
#'
#' @param ExpressionSet a standard PhyloExpressionSet or
#' DivergenceExpressionSet object.
#' @param permutations a numeric value specifying the number of permutations
#' that shall be performed. Default is \code{permutations} = \code{1000}.
#' @param parallel a boolean value specifying whether multicore (parallel)
#' processing shall be performed on a multicore machine.
#' @return a numeric \code{N x S} matrix containing Entropy
#' values for \code{S} developmental stages after \code{N} permutations.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{Entropy}}
#' @examples
#' 
#' # read example phylotranscriptomics data
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#' 
#' # perform N = 1000 sampling computations sequencially
#' sMatrix <- Sampled.Entropy.Matrix( ExpressionSet = PhyloExpressionSetExample, 
#'                                    permutations  = 1000 )
#' 
#' \dontrun{
#' 
#' # perform N = 1000 sampling computations on a multicore machine to speed up computations
#' psMatrix <- Sampled.Entropy.Matrix( ExpressionSet = PhyloExpressionSetExample, 
#'                                     permutations  = 1000,
#'                                     parallel      = TRUE )
#'}
#'
#' @import foreach
#' @export

Sampled.Entropy.Matrix <- function(ExpressionSet,permutations = 1000, parallel = FALSE){
  
    is.ExpressionSet(ExpressionSet)
    nCols <- dim(ExpressionSet)[2]
    nPS <- as.numeric(max(ExpressionSet[ , 1]))
    Entropy.Matrix <- matrix(NA_real_,permutations,(nCols-2))
    ExpressionMatrixConst <- as.matrix(ExpressionSet[ , 3:nCols]) 
    
    if(parallel){
      cores <- parallel::makeForkCluster(parallel::detectCores())
      doParallel::registerDoParallel(cores)
      
      ### Perform the sampling process in parallel
      foreach::foreach(i             = 1:permutations,
                      .combine       = "rbind",
                      .errorhandling = "stop") %dopar% {
                              
                              Entropy.Matrix[i ,] <- Phylotranscriptomics_Entropy(ExpressionMatrixConst,
                                                                                  as.vector(sample(ExpressionSet[ , 1])),
                                                                                  nPS)
                                }
      
      
      parallel::stopCluster(cores)
    
    }
    
    if(!parallel){
      ### perform computations sequencially 
      
      for(i in 1:permutations){
        
        Entropy.Matrix[i , ] <- Phylotranscriptomics_Entropy(ExpressionMatrixConst,as.vector(sample(ExpressionSet[ , 1])), nPS)
         
      }
         
      rownames(Entropy.Matrix) <- paste0("p",1:permutations)
      colnames(Entropy.Matrix) <- names(ExpressionSet)[3:nCols]
    
      return(Entropy.Matrix)
    
   }

}

