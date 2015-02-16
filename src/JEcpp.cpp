#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double JEcpp(NumericVector JointProbabilities){
  int len = JointProbabilities.size();
  double JointEntropy = 0;


     for(int i = 0; i < len; i++){
        if(JointProbabilities[i] > 0){

           JointEntropy += (JointProbabilities[i] * (log(JointProbabilities[i])/log(2)));

        } else{

          JointEntropy += 0;

        }
     }

  return(-JointEntropy);
}



