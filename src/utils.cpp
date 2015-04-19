#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
using namespace std;

//' @export
// [[Rcpp::export]]
double vectorSum(NumericVector x) {
        //if(is_na(x)){
        //       Rcpp::stop("Your input vector includes NA values...");
        // }
        
        // http://gallery.rcpp.org/articles/parallel-vector-sum/
        return std::accumulate(x.begin(), x.end(), 0.0);
}





