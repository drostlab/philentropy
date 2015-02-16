#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double MIcpp(NumericVector X, NumericVector Y, NumericVector XY){

  double MutualInformation = 0;
  // Using the identity: I(X,Y) = H(X) + H(Y) -  H(X,Y)
  MutualInformation = (Ecpp(X) + Ecpp(Y)) - JEcpp(XY);
  return(MutualInformation);
}


