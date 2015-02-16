#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
 NumericVector TAIcpp(NumericMatrix ExpressionSet, NumericVector Phylostratum){
    int nCols = ExpressionSet.ncol();
    int nRows = ExpressionSet.nrow();
    NumericVector results(nCols);
    for(int stage = 0; stage < nCols; stage++) {
      double numerator = 0, divisor = 0;
        for(int gene = 0; gene < nRows; gene++) {
		      numerator+= (double) Phylostratum[gene] * ExpressionSet(gene, stage);
		      divisor  += ExpressionSet(gene, stage);
	       }
	       results[stage] = numerator/divisor;
    }
    return results;
  }

