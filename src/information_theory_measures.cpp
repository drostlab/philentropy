#include <Rcpp.h>
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
                    Rcpp::stop("Please choose from units: log, log2, or log10.")
            }
    } else{
      Entropy += 0.0;
    }
  }

  return -Entropy;
}

//' @export
// [[Rcpp::export]]
double JEcpp(const NumericVector& JointProbabilities){
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
double CEcpp(const NumericVector& JointProbabilities,const NumericVector& Probabilities, const Rcpp::String unit){

  double ConditionalEntropy = 0.0;
  // Using the chain rule: H(X | Y) = H(X,Y) - H(Y)
  // Note: it is important that the Probabilities vector corresponds to Y
  ConditionalEntropy = JEcpp(JointProbabilities, unit) - Ecpp(Probabilities, unit);

  return(ConditionalEntropy);
}


//' @export
// [[Rcpp::export]]
double CrossEntropy(const NumericVector& P, const NumericVector& Q, const bool& testNA) {
     
     double log_2 = log(2);
     // Cross-Entropy = Kullback-Leibler Divergence
     double CE = 0.0;
     int Psize = P.size();
     int Qsize = Q.size();
     
     if(testNA){
               if(any(is_na(P)) | any(is_na(Q))){
                       Rcpp::stop("Your input vector stores NA values...");
                } 
     }
      
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
double MIcpp(const NumericVector& X, const NumericVector& Y, const NumericVector& XY){

  double MutualInformation = 0;
  // Using the identity: I(X,Y) = H(X) + H(Y) -  H(X,Y)
  MutualInformation = (Ecpp(X) + Ecpp(Y)) - JEcpp(XY);
  return(MutualInformation);
}


//' @export
// [[Rcpp::export]]
double JensonShannonDivergenceCpp(const NumericVector& P, const NumericVector& Q, const bool& testNA){
  int Psize = P.size();
  int Qsize = Q.size();
  double jsd = -1;

  if(Psize == Qsize){
    NumericVector R(Psize);

    for(int i=0; i < Psize; i++){
      R[i] = (P[i] + Q[i])/2.0;
    }

    jsd = 0.5 * (CrossEntropy(P,R,testNA) + CrossEntropy(Q,R,testNA));

  }

  return jsd;
}


