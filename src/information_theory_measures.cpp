#include <Rcpp.h>
#include "distances.h"

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace std;

//' @export
// [[Rcpp::export]]
double Ecpp(const NumericVector& P, const Rcpp::String unit){
  int len = P.size();
  double Entropy = 0.0;

  for (int i = 0; i < len; i++){
     if (P[i] > 0){
            if (unit == "log"){
                    Entropy += P[i] * log(P[i]);
            }
            
            else if (unit == "log2"){
                    Entropy += P[i] * custom_log2(P[i]);
            }
            
            else if (unit == "log10"){
                    Entropy += P[i] * custom_log10(P[i]);
                    
            } else {
                    Rcpp::stop("Please choose from units: log, log2, or log10.");
            }
    } else{
      Entropy += 0.0;
    }
  }

  return -Entropy;
}

//' @export
// [[Rcpp::export]]
double JEcpp(const NumericVector& JointProbabilities, const Rcpp::String unit){
  int len = JointProbabilities.size();
  double JointEntropy = 0.0;

     for (int i = 0; i < len; i++){
        if (JointProbabilities[i] > 0){
           
           if (unit == "log"){
                   JointEntropy += JointProbabilities[i] * log(JointProbabilities[i]);
           }
           
           else if (unit == "log2"){
                   JointEntropy += JointProbabilities[i] * custom_log2(JointProbabilities[i]);
           }
           
           else if (unit == "log10"){
                   JointEntropy += JointProbabilities[i] * custom_log10(JointProbabilities[i]);
                   
           } else {
                   Rcpp::stop("Please choose from units: log, log2, or log10.");
           }
        } else{

          JointEntropy += 0.0;

        }
     }

  return(-JointEntropy);
}


//' @export
// [[Rcpp::export]]
double CEcpp(const NumericVector& JointProbabilities,const NumericVector& Probabilities, const Rcpp::String unit){

  double ConditionalEntropy = 0.0;
  // Using the chain rule: H(X | Y) = H(X,Y) - H(Y)
  // Note: it is important that the Probabilities vector corresponds to Y
  ConditionalEntropy = JEcpp(JointProbabilities, unit) - Ecpp(Probabilities, unit);

  return(ConditionalEntropy);
}


//' @export
// [[Rcpp::export]]
double MIcpp(const NumericVector& X, const NumericVector& Y, const NumericVector& XY, const Rcpp::String unit){

  double MutualInformation = 0.0;
  // Using the identity: I(X,Y) = H(X) + H(Y) -  H(X,Y)
  MutualInformation = (Ecpp(X, unit) + Ecpp(Y, unit)) - JEcpp(XY, unit);
  
  return(MutualInformation);
}

