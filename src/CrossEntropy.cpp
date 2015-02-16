#include <Rcpp.h>
using namespace Rcpp;


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
                CE += (P[i] * (log((double) P[i]/Q[i])/log_2));
            } else{
              
                CE += 0;
         }
       }
        return CE;
     } else{       
        return -1;
    }
}