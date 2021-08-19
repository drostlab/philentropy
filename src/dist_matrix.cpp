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

//' @title Distances and Similarities between Two Probability Density Functions
//' @description This functions computes the distance/dissimilarity between two probability density functions.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param dist_fun a character string indicating whether the distance measure that should be computed.
//' @param p power of the Minkowski distance.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @param unit type of \code{log} function. Option are 
//' \itemize{
//' \item \code{unit = "log"}
//' \item \code{unit = "log2"}
//' \item \code{unit = "log10"}   
//' }
//' @param epsilon epsilon a small value to address cases in the distance computation where division by zero occurs. In
//' these cases, x / 0 or 0 / 0 will be replaced by \code{epsilon}. The default is \code{epsilon = 0.00001}.
//' However, we recommend to choose a custom \code{epsilon} value depending on the size of the input vectors,
//' the expected similarity between compared probability density functions and 
//' whether or not many 0 values are present within the compared vectors.
//' As a rough rule of thumb we suggest that when dealing with very large 
//' input vectors which are very similar and contain many \code{0} values,
//' the \code{epsilon} value should be set even smaller (e.g. \code{epsilon = 0.000000001}),
//' whereas when vector sizes are small or distributions very divergent then
//' higher \code{epsilon} values may also be appropriate (e.g. \code{epsilon = 0.01}).
//' Addressing this \code{epsilon} issue is important to avoid cases where distance metrics
//' return negative values which are not defined and only occur due to the
//' technical issues of computing x / 0 or 0 / 0 cases.
//' @return A single distance value
//' @export
// [[Rcpp::export]]
double dist_one_one(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, const Rcpp::String& dist_fun, const double& p = NA_REAL, const bool& testNA = true, const Rcpp::String& unit = "log", const double& epsilon = 0.00001){
        double dist_value;
        if (dist_fun == "euclidean"){
                dist_value = euclidean(P, Q, testNA);
        } else if (dist_fun == "manhattan") {
                dist_value = manhattan(P, Q, testNA);
        } else if (dist_fun == "minkowski") {
                dist_value = minkowski(P, Q, p, testNA);
        } else if (dist_fun == "chebyshev") {
                dist_value = chebyshev(P, Q, testNA);
        } else if (dist_fun == "sorensen") {
                dist_value = sorensen(P, Q, testNA);
        } else if (dist_fun == "gower") {
                dist_value = gower(P, Q, testNA);
        } else if (dist_fun == "soergel") {
                dist_value = soergel(P, Q, testNA);
        } else if (dist_fun == "kulczynski_d") {
                dist_value = kulczynski_d(P, Q, testNA, epsilon);
        } else if (dist_fun == "canberra") {
                dist_value = canberra(P, Q, testNA);
        } else if (dist_fun == "lorentzian") {
                dist_value = lorentzian(P, Q, testNA, unit);
        } else if (dist_fun == "intersection") {
                dist_value = intersection_dist(P, Q, testNA);
        } else if (dist_fun == "non-intersection") {
                dist_value = 1.0 - intersection_dist(P, Q, testNA);
        } else if (dist_fun == "wavehedges") {
                dist_value = wave_hedges(P, Q, testNA);
        } else if (dist_fun == "czekanowski") {
                dist_value = czekanowski(P, Q, testNA);
        } else if (dist_fun == "motyka") {
                dist_value = motyka(P, Q, testNA);
        } else if (dist_fun == "kulczynski_s") {
                dist_value = 1.0 / kulczynski_d(P, Q, testNA, epsilon);
        } else if (dist_fun == "tanimoto") {
                dist_value = tanimoto(P, Q, testNA);
        } else if (dist_fun == "ruzicka") {
                dist_value = ruzicka(P, Q, testNA);
        } else if (dist_fun == "inner_product") {
                dist_value = inner_product(P, Q, testNA);
        } else if (dist_fun == "harmonic_mean") {
                dist_value = harmonic_mean_dist(P, Q, testNA);
        } else if (dist_fun == "cosine") {
                dist_value = cosine_dist(P, Q, testNA);
        } else if (dist_fun == "hassebrook") {
                dist_value = kumar_hassebrook(P, Q, testNA);
        } else if (dist_fun == "jaccard") {
                dist_value = jaccard(P, Q, testNA);
        } else if (dist_fun == "dice") {
                dist_value = dice_dist(P, Q, testNA);
        } else if (dist_fun == "fidelity") {
                dist_value = fidelity(P, Q, testNA);
        } else if (dist_fun == "bhattacharyya") {
                dist_value = bhattacharyya(P, Q, testNA, unit, epsilon);
        } else if (dist_fun == "hellinger") {
                dist_value = hellinger(P, Q, testNA);
        } else if (dist_fun == "matusita") {
                dist_value = matusita(P, Q, testNA);
        } else if (dist_fun == "squared_chord") {
                dist_value = squared_chord(P, Q, testNA);
        } else if (dist_fun == "squared_euclidean") {
                dist_value = squared_euclidean(P, Q, testNA);
        } else if (dist_fun == "pearson") {
                dist_value = pearson_chi_sq(P, Q, testNA, epsilon);
        } else if (dist_fun == "neyman") {
                dist_value = neyman_chi_sq(P, Q, testNA, epsilon);
        } else if (dist_fun == "squared_chi") {
                dist_value = squared_chi_sq(P, Q, testNA);
        } else if (dist_fun == "prob_symm") {
                dist_value = prob_symm_chi_sq(P, Q, testNA);
        } else if (dist_fun == "divergence") {
                dist_value = divergence_sq(P, Q, testNA);
        } else if (dist_fun == "clark") {
                dist_value = clark_sq(P, Q, testNA);
        } else if (dist_fun == "additive_symm") {
                dist_value = additive_symm_chi_sq(P, Q, testNA);
        } else if (dist_fun == "kullback-leibler") {
                dist_value = kullback_leibler_distance(P, Q, testNA, unit, epsilon);
        } else if (dist_fun == "jeffreys") {
                dist_value = jeffreys(P, Q, testNA, unit, epsilon);
        } else if (dist_fun == "k_divergence") {
                dist_value = k_divergence(P, Q, testNA, unit);
        } else if (dist_fun == "topsoe") {
                dist_value = topsoe(P, Q, testNA, unit);
        } else if (dist_fun == "jensen_shannon"){
                dist_value = jensen_shannon(P, Q, testNA, unit); 
        } else if (dist_fun == "jensen_difference") {
                dist_value = jensen_difference(P, Q, testNA, unit);
        } else if (dist_fun == "taneja") {
                dist_value = taneja(P, Q, testNA, unit, epsilon);
        } else if (dist_fun == "kumar-johnson") {
                dist_value = kumar_johnson(P, Q, testNA, epsilon);
        } else if (dist_fun == "avg") {
                dist_value = avg(P, Q, testNA);
        } else {
                Rcpp::stop("Specified dist_fun is not implemented. Please consult getDistMethods().");
        }
        return dist_value;
}

//' @title Distances and Similarities between One and Many Probability Density Functions
//' @description This functions computes the distance/dissimilarity between one probability density functions and a set of probability density functions.
//' @param P a numeric vector storing the first distribution.
//' @param dists a numeric matrix storing distributions in its rows.
//' @param dist_fun a character string indicating whether the distance measure that should be computed.
//' @param p power of the Minkowski distance.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @param unit type of \code{log} function. Option are 
//' \itemize{
//' \item \code{unit = "log"}
//' \item \code{unit = "log2"}
//' \item \code{unit = "log10"}   
//' }
//' @param epsilon epsilon a small value to address cases in the distance computation where division by zero occurs. In
//' these cases, x / 0 or 0 / 0 will be replaced by \code{epsilon}. The default is \code{epsilon = 0.00001}.
//' However, we recommend to choose a custom \code{epsilon} value depending on the size of the input vectors,
//' the expected similarity between compared probability density functions and 
//' whether or not many 0 values are present within the compared vectors.
//' As a rough rule of thumb we suggest that when dealing with very large 
//' input vectors which are very similar and contain many \code{0} values,
//' the \code{epsilon} value should be set even smaller (e.g. \code{epsilon = 0.000000001}),
//' whereas when vector sizes are small or distributions very divergent then
//' higher \code{epsilon} values may also be appropriate (e.g. \code{epsilon = 0.01}).
//' Addressing this \code{epsilon} issue is important to avoid cases where distance metrics
//' return negative values which are not defined and only occur due to the
//' technical issues of computing x / 0 or 0 / 0 cases.
//' @return A vector of distance values
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dist_one_many(const Rcpp::NumericVector& P, Rcpp::NumericMatrix dists, Rcpp::String dist_fun, double p, bool testNA, Rcpp::String unit, double epsilon){
        
        int nrows = dists.nrow();
        Rcpp::NumericVector dist_values(nrows);
        
        for (int i = 0; i < nrows; i++){
                dist_values[i] = dist_one_one(P, dists(i, Rcpp::_), dist_fun, p, testNA, unit);
        }
        return dist_values;
}

//' @title Distances and Similarities between Many Probability Density Functions
//' @description This functions computes the distance/dissimilarity between two sets of probability density functions.
//' @param dists1 a numeric matrix storing distributions in its rows.
//' @param dists2 a numeric matrix storing distributions in its rows.
//' @param dist_fun a character string indicating whether the distance measure that should be computed.
//' @param p power of the Minkowski distance.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @param unit type of \code{log} function. Option are 
//' \itemize{
//' \item \code{unit = "log"}
//' \item \code{unit = "log2"}
//' \item \code{unit = "log10"}   
//' }
//' @param epsilon epsilon a small value to address cases in the distance computation where division by zero occurs. In
//' these cases, x / 0 or 0 / 0 will be replaced by \code{epsilon}. The default is \code{epsilon = 0.00001}.
//' However, we recommend to choose a custom \code{epsilon} value depending on the size of the input vectors,
//' the expected similarity between compared probability density functions and 
//' whether or not many 0 values are present within the compared vectors.
//' As a rough rule of thumb we suggest that when dealing with very large 
//' input vectors which are very similar and contain many \code{0} values,
//' the \code{epsilon} value should be set even smaller (e.g. \code{epsilon = 0.000000001}),
//' whereas when vector sizes are small or distributions very divergent then
//' higher \code{epsilon} values may also be appropriate (e.g. \code{epsilon = 0.01}).
//' Addressing this \code{epsilon} issue is important to avoid cases where distance metrics
//' return negative values which are not defined and only occur due to the
//' technical issues of computing x / 0 or 0 / 0 cases.
//' @return A matrix of distance values
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix dist_many_many(Rcpp::NumericMatrix dists1, Rcpp::NumericMatrix dists2, Rcpp::String dist_fun, double p, bool testNA, Rcpp::String unit, double epsilon){
        int nrows1 = dists1.nrow();
        int nrows2 = dists2.nrow();
        double dist_value = 0.0;
        
        Rcpp::NumericMatrix dist_matrix(nrows1,nrows2);
        // std::fill(dist_matrix.begin(), dist_matrix.end(), Rcpp::NumericVector::get_na());
        
        for (int i = 0; i < nrows1; i++){
                for (int j = 0; j < nrows2; j++){
                        // if(Rcpp::NumericVector::is_na(dist_matrix(i,j))){
                        dist_value = dist_one_one(dists1(i, Rcpp::_), dists2(j, Rcpp::_), dist_fun, p, testNA, unit);
                        dist_matrix(i,j) = dist_value;
                        // }
                }
        }
        return dist_matrix;
}
