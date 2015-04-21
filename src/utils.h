#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
using namespace std;

#ifndef philentropy_UTILS_H
#define philentropy_UTILS_H philentropy_UTILS_H

//[[Rcpp::export]]
Rcpp::NumericMatrix as_matrix(Rcpp::DataFrame x) {
// taken from: http://stackoverflow.com/questions/24352208/best-way-to-convert-dataframe-to-matrix-in-rcpp?rq=1        
  int nRows=x.nrows();  
  Rcpp::NumericMatrix y(nRows,x.size());
  for (int i=0; i<x.size();i++) {
       y(Rcpp::_,i) = Rcpp::NumericVector(x[i]);
  }  
  return y;
}

//' @export
// [[Rcpp::export]]
double vectorSum(NumericVector x) {
        //if(is_na(x)){
        //       Rcpp::stop("Your input vector includes NA values...");
        // }
        
        // http://gallery.rcpp.org/articles/parallel-vector-sum/
        return std::accumulate(x.begin(), x.end(), 0.0);
}

//' @export
// [[Rcpp::export]]
SEXP sum2( SEXP x_ ){
   NumericVector x(x_) ;
   double res = sum( x ) ;
   return wrap( res ) ;
}


 double sum3(NumericVector x) {
    double total = 0;
    
    NumericVector::iterator it;
    for(it = x.begin(); it != x.end(); ++it) {
      total += *it;
    }
    return total;
  }

#endif // philentropy_UTILS_H
