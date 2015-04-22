

#ifndef philentropy_UTILS_H
#define philentropy_UTILS_H philentropy_UTILS_H


#include <Rcpp.h>
#include <algorithm>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::interfaces(r, cpp)]]


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
Rcpp::DataFrame as_data_frame(Rcpp::NumericMatrix mat){
 
 //int nRows=mat.nrow();  
  Rcpp::DataFrame y;
  for (int i=0; i<mat.ncol();i++) {
       y[i] = mat(Rcpp::_,i);
  }  
  return y;
}


//' @export
// [[Rcpp::export]]
double vectorSum(Rcpp::NumericVector x) {
        //if(is_na(x)){
        //       Rcpp::stop("Your input vector includes NA values...");
        // }
        
        // http://gallery.rcpp.org/articles/parallel-vector-sum/
        return std::accumulate(x.begin(), x.end(), 0.0);
}

//' @export
// [[Rcpp::export]]
SEXP sum2( SEXP x_ ){
   Rcpp::NumericVector x(x_) ;
   double res = sum( x ) ;
   return Rcpp::wrap( res ) ;
}


 double sum3(Rcpp::NumericVector x) {
    double total = 0;
    
    Rcpp::NumericVector::iterator it;
    for(it = x.begin(); it != x.end(); ++it) {
      total += *it;
    }
    return total;
  }

#endif // philentropy_UTILS_H
