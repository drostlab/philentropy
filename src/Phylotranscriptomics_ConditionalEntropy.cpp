#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Phylotranscriptomics_ConditionalEntropy(NumericMatrix ExpressionSet, NumericVector X, NumericVector Y, int n_strata, int nY){
  // Conditional Entropy: H(X|Y) = H(Y,X) - H(Y)

  NumericVector H_yx = Phylotranscriptomics_JointEntropy(ExpressionSet,Y,X,nY,n_strata);
  NumericVector H_y = Phylotranscriptomics_Entropy(ExpressionSet,X,n_strata);
  int nStage = H_yx.size();
  NumericVector H_x_given_y(nStage);

  for(int stage = 0; stage < nStage; stage++){
      H_x_given_y[stage] = H_yx[stage] - H_y[stage];
    }

  return H_x_given_y;
}


