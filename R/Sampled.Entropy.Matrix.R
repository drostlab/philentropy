#' @title Perform Sampling for Entropy
#' 
#' @description
#'This function receives a standard PhyloExpressionSet or
#' DivergenceExpressionSet, samples the PS or DS column of the corresponding
#' PhyloExpressionSet or DivergenceExpressionSet and furthermore, computes the
#' phylotranscriptomics Entropy profile of this sampled PhyloExpressionSet or
#' DivergenceExpressionSet. This procedure is being done N times, with N =
#' number of permutations.
#' 
#' The resulting sample matrix consists of N x S elements, where S denotes the
#' total number of developmental stages stored in the corresponding
#' PhyloExpressionSet or DivergenceExpressionSet. This N x S sample matrix can
#' be used to perform further test statistics (sampling test) of the
#' statistical significance of any phylotranscriptomics Entropy pattern.
#' 
#' To speed up computations the sampling process can be paralallized (on a
#' multicore machine).
#' 
#' Internally this functions calls the C++ function
#' \code{\link{Phylotranscriptomics_Entropy}} to speed up the sampling process
#' and furthermore (optionally) uses muticore computing.
#' 
#' @param ExpressionSet a standard PhyloExpressionSet or
#' DivergenceExpressionSet object.
#' @param permutations a numeric value specifying the number of permutations
#' that shall be performed.
#' @param parallel a boolean value specifying whether multicore (parallel)
#' computations shall be performed on a multicore machine. The parallelization
#' only works on a single multicore machine and not on a cluster.
#' @return a numeric N x S matrix containing phylotranscriptomics Entropy
#' values for S developmental stages after N samplings.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{Entropy}}
#' @examples
#' 
#' # read standard phylotranscriptomics data
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#' 
#' # perform N=1000 sampling computations sequencially
#' sMatrix <- Sampled.Entropy.Matrix(PhyloExpressionSetExample,nPermutations=1000)
#' 
#' # perform N=1000 sampling computations on a multicore machine to speed up computations
#' psMatrix <- Sampled.Entropy.Matrix(PhyloExpressionSetExample,nPermutations=1000,parallel=T)
#' @export

Sampled.Entropy.Matrix <- function(ExpressionSet,permutations = 1000, parallel = FALSE){
  
    is.ExpressionSet(ExpressionSet)
    nCols <- dim(ExpressionSet)[2]
    nPS <- as.numeric(max(ExpressionSet[,1]))
    Entropy.Matrix <- matrix(0,nPermutations,(nCols-2))
    ExpressionMatrixConst <- as.matrix(ExpressionSet[,3:nCols]) 
    
    if(parallel){
      ### Parallellizing the sampling process using the 'doMC' and 'parallel' package
      ### register all given cores for parallelization
      ### detectCores(all.tests = TRUE, logical = FALSE) returns the number of cores available on a multi-core machine
      cores=makeForkCluster(detectCores(all.tests = FALSE, logical = FALSE))
      registerDoParallel(cores)
      
      ### Perform the sampling process in parallel
      Entropy.Matrix <- as.matrix(foreach(i = 1:nPermutations,.combine="rbind") %dopar% {Phylotranscriptomics_Entropy(ExpressionMatrixConst,as.vector(sample(ExpressionSet[,1])),nPS)})
      
      ### close the cluster connection
      ### The is important to be able to re-run the function N times
      ### without getting cluster connection problems
      stopCluster(cores)
    }
    
    if(!parallel){
      ### perform computations sequencially 
      
      if(nPermutations > 100){
        ### initializing the progress bar
        progressBar <- txtProgressBar(min=1,max=nPermutations,style=3)
      }
      for(i in 1:nPermutations){
        
        Entropy.Matrix[i,] <- Phylotranscriptomics_Entropy(ExpressionMatrixConst,as.vector(sample(ExpressionSet[,1])),nPS)
        
        if(nPermutations > 100){
          ### printing out the progress
          setTxtProgressBar(progressBar,i)
        }
      }
      print("\n")
    }
    
    rownames(Entropy.Matrix) <- 1:nPermutations
    colnames(Entropy.Matrix) <- names(ExpressionSet)[3:nCols]
    
    return(Entropy.Matrix)
    
}

