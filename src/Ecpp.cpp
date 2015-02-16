#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double Ecpp(NumericVector Probabilities){
  int len = Probabilities.size();
  double Entropy = 0;

  for(int i = 0; i < len; i++){
    if(Probabilities[i] > 0){
      Entropy += (Probabilities[i] * (log(Probabilities[i])/log(2)));
    }

    else{
      Entropy += 0;
    }
  }

  return(-Entropy);
}

