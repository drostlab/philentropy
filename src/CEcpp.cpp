#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double CEcpp(NumericVector JointProbabilities,NumericVector Probabilities){

  double ConditionalEntropy = 0;
  // Using the chain rule: H(X | Y) = H(X,Y) - H(Y)
  // Note: it is important that the Probabilities vector corresponds to Y
  ConditionalEntropy = JEcpp(JointProbabilities) - Ecpp(Probabilities);

  return(ConditionalEntropy);
}

