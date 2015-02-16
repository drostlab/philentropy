#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Phylotranscriptomics_JointEntropy(NumericMatrix ExpressionSet, NumericVector X, NumericVector Y, int n_strata, int nY) {
   double log_2 = log(2);
   int nCols = ExpressionSet.ncol();
   NumericMatrix P_xy = Phylotranscriptomics_JointProbability(ExpressionSet,X,Y,n_strata,nY);
   int nProbabilities = P_xy.nrow();
   NumericVector joint_entropy(nCols);

   for(int stage = 0; stage < nCols; stage++){
      for(int it = 0; it < nProbabilities; it++){
          if(P_xy(it,stage) > 0){
             joint_entropy[stage] += (P_xy(it,stage) * ((double) log(P_xy(it,stage))/log_2));
          }
          else{
             joint_entropy[stage] += 0;
          }
      }
      joint_entropy[stage] *= -1;
    }

  return joint_entropy;
}


