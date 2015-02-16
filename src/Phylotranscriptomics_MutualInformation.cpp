#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Phylotranscriptomics_MutualInformation(NumericMatrix ExpressionSet, NumericVector X, NumericVector Y, int n_strata, int nY) {
   // Mutual Information : I(X,Y) = H(X) + H(Y) - H(X,Y)
   NumericVector H_x = Phylotranscriptomics_Entropy(ExpressionSet,X,n_strata);
   NumericVector H_y = Phylotranscriptomics_Entropy(ExpressionSet,Y,nY);
   NumericVector H_xy = Phylotranscriptomics_JointEntropy(ExpressionSet,X,Y,n_strata,nY);
   int nStages = H_x.size();
   NumericVector mutual_information(nStages);

   for(int stage= 0; stage < nStages; stage++){
       mutual_information[stage] = (H_x[stage] + H_y[stage]) - H_xy[stage];
   }

  return mutual_information;

}

