#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix Phylotranscriptomics_JointProbability(NumericMatrix ExpressionSet, NumericVector X, NumericVector Y, int n_strata, int nY){

   int nCols = ExpressionSet.ncol();
   int nRows = ExpressionSet.nrow();
   int iterator = 0;
   NumericMatrix resultMatrix(n_strata*nY,nCols);

   for(int stage = 0; stage < nCols; stage++){
     for(int x_it = 1; x_it <= n_strata; x_it++){
       for(int y_it = 1; y_it <= nY; y_it++){

         double numerator = 0, divisor = 0;

         for(int gene = 0; gene < nRows; gene++){

            if((X[gene] == x_it) && (Y[gene] == y_it)){

               numerator += ExpressionSet(gene,stage);

            }

            divisor  += ExpressionSet(gene,stage);
         }
         iterator = (y_it-1) + ((x_it-1)*nY);
         resultMatrix(iterator,stage) = (double) numerator/divisor;
      }
   }
}
return(resultMatrix);
}


