#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Phylotranscriptomics_Entropy(NumericMatrix ExpressionSet, NumericVector X, int n_strata) {
     double log_2 = log(2);
     int nCols = ExpressionSet.ncol();
     NumericVector entropy(nCols);
     NumericMatrix P_x = Phylotranscriptomics_Probability(ExpressionSet, X, n_strata);

     for(int stage = 0; stage < nCols; stage++){
          for(int x_it = 0; x_it < n_strata; x_it++){
               if(P_x(x_it,stage) > 0){
                   entropy[stage] += (P_x(x_it,stage) * ((double) log(P_x(x_it,stage))/log_2));
               } else{

                   entropy[stage] += 0;

               }
          }
      entropy[stage] *= -1;

     }

     return entropy;
}