// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include "utils.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix as_matrix(Rcpp::DataFrame x) {
// taken from: http://stackoverflow.com/questions/24352208/best-way-to-convert-dataframe-to-matrix-in-rcpp?rq=1
  int nRows = x.nrows();
  Rcpp::NumericMatrix y(nRows, x.size());
  for (int i = 0; i < x.size(); i++) {
       y(Rcpp::_, i) = Rcpp::NumericVector(x[i]);
  }
  return y;
}

// @export
// [[Rcpp::export]]
Rcpp::DataFrame as_data_frame(Rcpp::NumericMatrix mat) {
  Rcpp::DataFrame y;
  for (int i = 0; i < mat.ncol(); i++) {
       y[i] = mat(Rcpp::_, i);
  }
  return y;
}

// @export
// [[Rcpp::export]]
SEXP sum_rcpp(SEXP vec) {
   Rcpp::NumericVector x(vec);
   double res = sum(x);
   return Rcpp::wrap(res);
}

// @export
// [[Rcpp::export]]
SEXP est_prob_empirical(SEXP CountVec) {
   Rcpp::NumericVector x(CountVec);
   double ProbMass = sum(x);
   Rcpp::NumericVector EmpiricalProb(x.size());

   EmpiricalProb = x / ProbMass;

   return Rcpp::wrap(EmpiricalProb);
}