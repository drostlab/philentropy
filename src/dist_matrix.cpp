// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "utils.h"
#include "distances.h"
#include "InformationTheory.h"
#include "correlation.h"

//using namespace Rcpp;
//using namespace std;


// @export
// [[Rcpp::export]]
Rcpp::NumericMatrix DistMatrixWithoutUnitDF(Rcpp::DataFrame distsDF, Rcpp::Function DistFunc, bool testNA){
// http://stackoverflow.com/questions/27391472/passing-r-function-as-parameter-to-rcpp-function
        Rcpp::NumericMatrix dists = as_matrix(distsDF);
        int nrow = dists.nrow();
        double dist_value = 0.0;
        Rcpp::NumericMatrix dist_matrix(nrow,nrow);
        // http://stackoverflow.com/questions/23748572/initializing-a-matrix-to-na-in-rcpp
        std::fill( dist_matrix.begin(), dist_matrix.end(), Rcpp::NumericVector::get_na() );
        
        for (int i = 0; i < nrow; i++){
                for (int j = 0; j < nrow; j++){
                        if(Rcpp::NumericVector::is_na(dist_matrix(i,j))){
                                dist_value = Rcpp::as<double>(DistFunc(dists(i,Rcpp::_),dists(j,Rcpp::_), testNA));
                                dist_matrix(i,j) = dist_value;
                                dist_matrix(j,i) = dist_value;
                        }
                }
        }
        
        return dist_matrix;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix DistMatrixMinkowskiMAT(Rcpp::NumericMatrix dists, double p, bool testNA){
        int ncols = dists.ncol();
        double dist_value = 0.0;
        Rcpp::NumericMatrix dist_matrix(ncols,ncols);
        // http://stackoverflow.com/questions/23748572/initializing-a-matrix-to-na-in-rcpp
        std::fill( dist_matrix.begin(), dist_matrix.end(), Rcpp::NumericVector::get_na() );
        
        for (int i = 0; i < ncols; i++){
                for (int j = 0; j < ncols; j++){
                        if(Rcpp::NumericVector::is_na(dist_matrix(i,j))){
                                dist_value = minkowski(dists(Rcpp::_, i),dists(Rcpp::_, j), p, testNA);
                                dist_matrix(i,j) = dist_value;
                                dist_matrix(j,i) = dist_value;
                        }
                }
        }
        
        return dist_matrix;
}

// @export
// [[Rcpp::export]]
Rcpp::NumericMatrix DistMatrixWithoutUnitMAT(Rcpp::NumericMatrix dists, Rcpp::Function DistFunc, bool testNA){
        // http://stackoverflow.com/questions/27391472/passing-r-function-as-parameter-to-rcpp-function
        
        int ncols = dists.ncol();
        double dist_value = 0.0;
        Rcpp::NumericMatrix dist_matrix(ncols,ncols);
        // http://stackoverflow.com/questions/23748572/initializing-a-matrix-to-na-in-rcpp
        std::fill( dist_matrix.begin(), dist_matrix.end(), Rcpp::NumericVector::get_na() );
        
        for (int i = 0; i < ncols; i++){
                for (int j = 0; j < ncols; j++){
                        if(Rcpp::NumericVector::is_na(dist_matrix(i,j))){
                                dist_value = Rcpp::as<double>(DistFunc(dists(Rcpp::_, i),dists(Rcpp::_, j), testNA));
                                dist_matrix(i,j) = dist_value;
                                dist_matrix(j,i) = dist_value;
                        }
                }
        }
        
        return dist_matrix;
}


// @export
// [[Rcpp::export]]
Rcpp::NumericMatrix DistMatrixWithUnitDF(Rcpp::DataFrame distsDF, Rcpp::Function DistFunc, bool testNA, Rcpp::String unit){
// http://stackoverflow.com/questions/27391472/passing-r-function-as-parameter-to-rcpp-function

        Rcpp::NumericMatrix dists = as_matrix(distsDF);
        int nrow = dists.nrow();
        double dist_value = 0.0;
        Rcpp::NumericMatrix dist_matrix(nrow,nrow);
        // http://stackoverflow.com/questions/23748572/initializing-a-matrix-to-na-in-rcpp
        std::fill( dist_matrix.begin(), dist_matrix.end(), Rcpp::NumericVector::get_na() );
        
        for (int i = 0; i < nrow; i++){
                for (int j = 0; j < nrow; j++){
                        if(Rcpp::NumericVector::is_na(dist_matrix(i,j))){
                                dist_value = Rcpp::as<double>(DistFunc(dists(i,Rcpp::_),dists(j,Rcpp::_), testNA, unit));
                                dist_matrix(i,j) = dist_value;
                                dist_matrix(j,i) = dist_value;
                        }
                }
        }
        
        return dist_matrix;
}


// @export
// [[Rcpp::export]]
Rcpp::NumericMatrix DistMatrixWithUnitMAT(Rcpp::NumericMatrix dists, Rcpp::Function DistFunc, bool testNA, Rcpp::String unit){
// http://stackoverflow.com/questions/27391472/passing-r-function-as-parameter-to-rcpp-function

        int ncols = dists.ncol();
        double dist_value = 0.0;
        Rcpp::NumericMatrix dist_matrix(ncols,ncols);
        // http://stackoverflow.com/questions/23748572/initializing-a-matrix-to-na-in-rcpp
        std::fill( dist_matrix.begin(), dist_matrix.end(), Rcpp::NumericVector::get_na() );
        
        for (int i = 0; i < ncols; i++){
                for (int j = 0; j < ncols; j++){
                        if(Rcpp::NumericVector::is_na(dist_matrix(i,j))){
                                dist_value = Rcpp::as<double>(DistFunc(dists(Rcpp::_, i),dists(Rcpp::_, j), testNA, unit));
                                dist_matrix(i,j) = dist_value;
                                dist_matrix(j,i) = dist_value;
                        }
                }
        } 
        
        return dist_matrix;
}

// @export
// [[Rcpp::export]]
double single_distance(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, const Rcpp::String& dist_fun, const bool& testNA, const Rcpp::String& unit){
        double dist_value;
        if (dist_fun == "euclidean"){
                dist_value = euclidean(P, Q, testNA);
        } else if (dist_fun == "jensen_shannon"){
                dist_value = jensen_shannon(P, Q, testNA, unit);
        }
        return dist_value;
}

// @export
// [[Rcpp::export]]
Rcpp::NumericVector dist_one_many(const Rcpp::NumericVector& P, Rcpp::NumericMatrix dists, Rcpp::String dist_fun, bool testNA, Rcpp::String unit){
        
        int nrows = dists.nrow();
        Rcpp::NumericVector dist_values(nrows);
        
        for (int i = 0; i < nrows; i++){
                dist_values[i] = single_distance(P, dists(i, Rcpp::_), dist_fun, testNA, unit);
        }
        return dist_values;
}

// @export
// [[Rcpp::export]]
Rcpp::NumericMatrix dist_many_many(Rcpp::NumericMatrix dists1, Rcpp::NumericMatrix dists2, Rcpp::String dist_fun, bool testNA, Rcpp::String unit){
        int nrows1 = dists1.nrow();
        int nrows2 = dists2.nrow();
        double dist_value = 0.0;
        
        Rcpp::NumericMatrix dist_matrix(nrows1,nrows2);
        // std::fill(dist_matrix.begin(), dist_matrix.end(), Rcpp::NumericVector::get_na());
        
        for (int i = 0; i < nrows1; i++){
                for (int j = 0; j < nrows2; j++){
                        // if(Rcpp::NumericVector::is_na(dist_matrix(i,j))){
                                dist_value = single_distance(dists1(i, Rcpp::_), dists2(j, Rcpp::_), dist_fun, testNA, unit);
                                dist_matrix(i,j) = dist_value;
                        // }
                }
        }
        return dist_matrix;
}
