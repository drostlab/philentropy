#include <Rcpp.h>
using namespace Rcpp;

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

