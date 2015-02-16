#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix Phylotranscriptomics_Probability(NumericMatrix ExpressionSet, NumericVector X, int n_strata) {

      int nCols = ExpressionSet.ncol();
      int nRows = ExpressionSet.nrow();
      NumericMatrix results(n_strata,nCols);
      
        for(int stage = 0; stage < nCols; stage++) {
           for(int ps = 1; ps <= n_strata; ps++){

            double numerator = 0, divisor = 0;

             for(int gene = 0; gene < nRows; gene++) {

                  if(X[gene] == ps){

                      numerator += ExpressionSet(gene, stage);

                  }

                  divisor  += ExpressionSet(gene, stage);
             }

             results(ps-1,stage) = (double) numerator/divisor;
     }
    }
    return results;
}


