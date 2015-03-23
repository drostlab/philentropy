#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace std;

//' @export
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

//' @export
// [[Rcpp::export]]
double JEcpp(NumericVector JointProbabilities){
  int len = JointProbabilities.size();
  double JointEntropy = 0;

     for(int i = 0; i < len; i++){
        if(JointProbabilities[i] > 0){

           JointEntropy += (JointProbabilities[i] * (log(JointProbabilities[i])/log(2)));

        } else{

          JointEntropy += 0.0;

        }
     }

  return(-JointEntropy);
}


//' @export
// [[Rcpp::export]]
double CEcpp(NumericVector JointProbabilities,NumericVector Probabilities){

  double ConditionalEntropy = 0;
  // Using the chain rule: H(X | Y) = H(X,Y) - H(Y)
  // Note: it is important that the Probabilities vector corresponds to Y
  ConditionalEntropy = JEcpp(JointProbabilities) - Ecpp(Probabilities);

  return(ConditionalEntropy);
}


//' @export
// [[Rcpp::export]]
double CrossEntropy(NumericVector P, NumericVector Q) {
     double log_2 = log(2);
     int Psize = P.size();
     int Qsize = Q.size();
     
     // Cross-Entropy = Kullback-Leibler Divergence
     double CE = 0.0;
     if(Psize == Qsize){
         for(int i = 0; i < Psize; i++){
            if((P[i] > 0) && (Q[i] > 0)){
                    
                CE += (P[i] * (log(P[i]/Q[i])/log_2));
                
            } else{
              
                CE += 0.0;
         }
       }
        return CE;
        
     } else{       
             
        return -1;
    }
}

//' @export
// [[Rcpp::export]]
double MIcpp(NumericVector X, NumericVector Y, NumericVector XY){

  double MutualInformation = 0;
  // Using the identity: I(X,Y) = H(X) + H(Y) -  H(X,Y)
  MutualInformation = (Ecpp(X) + Ecpp(Y)) - JEcpp(XY);
  return(MutualInformation);
}


//' @export
// [[Rcpp::export]]
double JensonShannonDivergenceCpp(NumericVector P,NumericVector Q){
  int Psize = P.size();
  int Qsize = Q.size();
  double jsd = -1;

  if(Psize == Qsize){
    NumericVector R(Psize);

    for(int i=0; i < Psize; i++){
      R[i] = (P[i] + Q[i])/2;
    }

    jsd = 0.5 * (CrossEntropy(P,R) + CrossEntropy(Q,R));

  }

  return jsd;
}


