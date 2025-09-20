// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include "utils.h"
#include "distances_internal.h"
#include "distances.h"
#include "InformationTheory.h"
#include "correlation.h"
#include <RcppParallel.h>

//using namespace Rcpp;
//using namespace std;

// Forward declarations for internal parallel functions
Rcpp::NumericMatrix DistMatrixMinkowskiMAT_internal(Rcpp::NumericMatrix dists, double p, bool testNA);
Rcpp::NumericMatrix DistMatrixWithoutUnitMAT_internal(Rcpp::NumericMatrix dists, std::string method, bool testNA, Rcpp::Nullable<double> p);
Rcpp::NumericMatrix DistMatrixWithUnitMAT_internal(Rcpp::NumericMatrix dists, std::string method, bool testNA, std::string unit);
Rcpp::NumericMatrix DistMatrixWithoutUnitDF_internal(Rcpp::DataFrame distsDF, std::string method, bool testNA, Rcpp::Nullable<double> p);
Rcpp::NumericMatrix DistMatrixWithUnitDF_internal(Rcpp::DataFrame distsDF, std::string method, bool testNA, std::string unit);
Rcpp::NumericMatrix dist_many_many_parallel(Rcpp::NumericMatrix& dists1, Rcpp::NumericMatrix& dists2, Rcpp::String method, double p, bool testNA, Rcpp::String unit, double epsilon);

// [[Rcpp::export]]
Rcpp::NumericMatrix DistMatrixWithoutUnitDF(Rcpp::DataFrame distsDF,
                                            std::string DistFunc,
                                            bool testNA,
                                            Rcpp::Nullable<double> p = R_NilValue) {
    return DistMatrixWithoutUnitDF_internal(distsDF, DistFunc, testNA, p);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix DistMatrixMinkowskiMAT(Rcpp::NumericMatrix dists, double p, bool testNA){
    return DistMatrixMinkowskiMAT_internal(dists, p, testNA);
}

struct MinkowskiWorker : public RcppParallel::Worker {
    RcppParallel::RMatrix<double> dists;
    RcppParallel::RMatrix<double> dist_matrix;
    double p;

    MinkowskiWorker(Rcpp::NumericMatrix& dists,
                    Rcpp::NumericMatrix& dist_matrix,
                       double p)
        : dists(dists), dist_matrix(dist_matrix), p(p) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            for (std::size_t j = i; j < (std::size_t)dists.ncol(); ++j) {
                auto col_i = dists.column(i);
                auto col_j = dists.column(j);
                double dist = minkowski_internal(col_i.begin(), col_i.end(), col_j.begin(), p);
                dist_matrix(i, j) = dist;
                dist_matrix(j, i) = dist;
            }
        }
    }
};

Rcpp::NumericMatrix DistMatrixMinkowskiMAT_internal(Rcpp::NumericMatrix dists,
                                                    double p, bool testNA) {
    int n = dists.ncol();
    Rcpp::NumericMatrix dist_matrix(n, n);

    MinkowskiWorker worker(dists, dist_matrix, p);
    RcppParallel::parallelFor(0, n, worker);

    return dist_matrix;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix DistMatrixWithoutUnitMAT(Rcpp::NumericMatrix dists,
                                             std::string DistFunc,
                                             bool testNA,
                                             Rcpp::Nullable<double> p = R_NilValue) {
    return DistMatrixWithoutUnitMAT_internal(dists, DistFunc, testNA, p);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix DistMatrixWithUnitDF(Rcpp::DataFrame distsDF, std::string DistFunc, bool testNA, std::string unit) {
    return DistMatrixWithUnitDF_internal(distsDF, DistFunc, testNA, unit);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix DistMatrixWithUnitMAT(Rcpp::NumericMatrix dists, std::string DistFunc, bool testNA, std::string unit) {
    return DistMatrixWithUnitMAT_internal(dists, DistFunc, testNA, unit);
}

//' @title Distances and Similarities between Two Probability Density Functions
//' @description This functions computes the distance/dissimilarity between two probability density functions.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param method a character string indicating whether the distance measure that should be computed.
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
//' @examples
//' P <- 1:10 / sum(1:10)
//' Q <- 20:29 / sum(20:29)
//' dist_one_one(P, Q, method = "euclidean", testNA = FALSE)
//' @export
// [[Rcpp::export]]
double dist_one_one(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, const Rcpp::String& method, Rcpp::Nullable<double> p = R_NilValue, const bool& testNA = true, const Rcpp::String& unit = "log", const double& epsilon = 0.00001){
        double dist_value;
        std::string unit_str(unit);

        if (method == "euclidean"){
                dist_value = euclidean_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "manhattan") {
                dist_value = minkowski_internal(P.begin(), P.end(), Q.begin(), 1.0);
        } else if (method == "minkowski") {
                if (p.isNotNull()) {
                    dist_value = minkowski_internal(P.begin(), P.end(), Q.begin(), Rcpp::as<double>(p));
                } else {
                    Rcpp::stop("Please specify p for the Minkowski distance.");
                }
        } else if (method == "chebyshev") {
                dist_value = chebyshev_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "sorensen") {
                dist_value = sorensen_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "gower") {
                dist_value = gower_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "soergel") {
                dist_value = soergel_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "kulczynski_d") {
                dist_value = kulczynski_d_internal(P.begin(), P.end(), Q.begin(), epsilon);
        } else if (method == "canberra") {
                dist_value = canberra_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "lorentzian") {
                dist_value = lorentzian_internal(P.begin(), P.end(), Q.begin(), unit_str);
        } else if (method == "intersection") {
                dist_value = intersection_dist_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "non-intersection") {
                dist_value = 1.0 - intersection_dist_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "wavehedges") {
                dist_value = wave_hedges_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "czekanowski") {
                dist_value = czekanowski_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "motyka") {
                dist_value = motyka_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "kulczynski_s") {
                dist_value = 1.0 / kulczynski_d_internal(P.begin(), P.end(), Q.begin(), epsilon);
        } else if (method == "tanimoto") {
                dist_value = tanimoto_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "ruzicka") {
                dist_value = ruzicka_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "inner_product") {
                dist_value = inner_product_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "harmonic_mean") {
                dist_value = harmonic_mean_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "cosine") {
                dist_value = cosine_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "hassebrook") {
                dist_value = hassebrook_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "jaccard") {
                dist_value = jaccard_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "dice") {
                dist_value = dice_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "fidelity") {
                dist_value = fidelity_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "bhattacharyya") {
                dist_value = bhattacharyya_internal(P.begin(), P.end(), Q.begin(), unit_str, epsilon);
        } else if (method == "hellinger") {
                dist_value = hellinger_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "matusita") {
                dist_value = matusita_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "squared_chord") {
                dist_value = squared_chord_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "squared_euclidean") {
                dist_value = squared_euclidean_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "pearson") {
                dist_value = pearson_internal(P.begin(), P.end(), Q.begin(), epsilon);
        } else if (method == "neyman") {
                dist_value = neyman_internal(P.begin(), P.end(), Q.begin(), epsilon);
        } else if (method == "squared_chi") {
                dist_value = squared_chi_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "prob_symm") {
                dist_value = prob_symm_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "divergence") {
                dist_value = divergence_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "clark") {
                dist_value = clark_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "additive_symm") {
                dist_value = additive_symm_internal(P.begin(), P.end(), Q.begin());
        } else if (method == "kullback-leibler") {
                dist_value = kullback_leibler_internal(P.begin(), P.end(), Q.begin(), unit_str, epsilon);
        } else if (method == "jeffreys") {
                dist_value = jeffreys_internal(P.begin(), P.end(), Q.begin(), unit_str, epsilon);
        } else if (method == "k_divergence") {
                dist_value = k_divergence_internal(P.begin(), P.end(), Q.begin(), unit_str);
        } else if (method == "topsoe") {
                dist_value = topsoe_internal(P.begin(), P.end(), Q.begin(), unit_str);
        } else if (method == "jensen-shannon"){
                dist_value = jensen_shannon_internal(P.begin(), P.end(), Q.begin(), unit_str); 
        } else if (method == "jensen_difference") {
                dist_value = jensen_difference_internal(P.begin(), P.end(), Q.begin(), unit_str);
        } else if (method == "taneja") {
                dist_value = taneja_internal(P.begin(), P.end(), Q.begin(), unit_str, epsilon);
        } else if (method == "kumar-johnson") {
                dist_value = kumar_johnson_internal(P.begin(), P.end(), Q.begin(), epsilon);
        } else if (method == "avg") {
                dist_value = avg_internal(P.begin(), P.end(), Q.begin());
        } else {
                Rcpp::stop("Specified method is not implemented. Please consult getDistMethods().");
        }
        return dist_value;
}

//' @title Distances and Similarities between One and Many Probability Density Functions
//' @description This functions computes the distance/dissimilarity between one probability density functions and a set of probability density functions.
//' @param P a numeric vector storing the first distribution.
//' @param dists a numeric matrix storing distributions in its rows.
//' @param method a character string indicating whether the distance measure that should be computed.
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
//' @examples
//' set.seed(2020-08-20)
//' P <- 1:10 / sum(1:10)
//' M <- t(replicate(100, sample(1:10, size = 10) / 55))
//' dist_one_many(P, M, method = "euclidean", testNA = FALSE)
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dist_one_many(const Rcpp::NumericVector& P, Rcpp::NumericMatrix dists, Rcpp::String method, Rcpp::Nullable<double> p = R_NilValue, bool testNA = true, Rcpp::String unit = "log", double epsilon = 0.00001){
        
        int nrows = dists.nrow();
        Rcpp::NumericVector dist_values(nrows);
        R_xlen_t vector_elements = nrows;
        
        
        for (R_xlen_t i = 0; i < vector_elements; i++){
                dist_values[i] = dist_one_one(P, dists(i, Rcpp::_), method, p, testNA, unit, epsilon);
        }
        return dist_values;
}

//' @title Distances and Similarities between Many Probability Density Functions
//' @description This functions computes the distance/dissimilarity between two sets of probability density functions.
//' @param dists1 a numeric matrix storing distributions in its rows.
//' @param dists2 a numeric matrix storing distributions in its rows.
//' @param method a character string indicating whether the distance measure that should be computed.
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
//' @examples 
//'   set.seed(2020-08-20)
//'   M1 <- t(replicate(10, sample(1:10, size = 10) / 55))
//'   M2 <- t(replicate(10, sample(1:10, size = 10) / 55))
//'   result <- dist_many_many(M1, M2, method = "euclidean", testNA = FALSE)
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix dist_many_many(Rcpp::NumericMatrix dists1, Rcpp::NumericMatrix dists2, Rcpp::String method, double p = NA_REAL, bool testNA = true, Rcpp::String unit = "log", double epsilon = 0.00001) {
    return dist_many_many_parallel(dists1, dists2, method, p, testNA, unit, epsilon);
}

struct ManyManyWorker : public RcppParallel::Worker {
    RcppParallel::RMatrix<double> dists1;
    RcppParallel::RMatrix<double> dists2;
    RcppParallel::RMatrix<double> dist_matrix;
    std::string method;
    double p;
    std::string unit;
    double epsilon;

    ManyManyWorker(Rcpp::NumericMatrix& dists1,
                      Rcpp::NumericMatrix& dists2,
                      Rcpp::NumericMatrix& dist_matrix,
                      std::string method,
                      double p,
                      std::string unit,
                      double epsilon)
        : dists1(dists1), dists2(dists2), dist_matrix(dist_matrix),
          method(method), p(p), unit(unit), epsilon(epsilon) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            for (std::size_t j = 0; j < (std::size_t)dists2.nrow(); ++j) {
                auto row_i = dists1.row(i);
                auto row_j = dists2.row(j);
                double dist = 0;
                if (method == "euclidean") {
                    dist = euclidean_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "manhattan") {
                    dist = minkowski_internal(row_i.begin(), row_i.end(), row_j.begin(), 1.0);
                } else if (method == "minkowski") {
                    dist = minkowski_internal(row_i.begin(), row_i.end(), row_j.begin(), p);
                } else if (method == "chebyshev") {
                    double max_dist = 0.0;
                    for(size_t k = 0; k < (size_t)row_i.length(); ++k) {
                        max_dist = std::max(max_dist, std::abs(row_i[k] - row_j[k]));
                    }
                    dist = max_dist;
                } else if (method == "sorensen") {
                    dist = sorensen_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "gower") {
                    dist = gower_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "soergel") {
                    dist = soergel_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "kulczynski_d") {
                    dist = kulczynski_d_internal(row_i.begin(), row_i.end(), row_j.begin(), epsilon);
                } else if (method == "canberra") {
                    dist = canberra_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "lorentzian") {
                    dist = lorentzian_internal(row_i.begin(), row_i.end(), row_j.begin(), unit);
                } else if (method == "intersection") {
                    dist = intersection_dist_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "wavehedges") {
                    dist = wave_hedges_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "czekanowski") {
                    dist = czekanowski_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "motyka") {
                    dist = motyka_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "tanimoto") {
                    dist = tanimoto_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "ruzicka") {
                    dist = ruzicka_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "inner_product") {
                    dist = inner_product_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "harmonic_mean") {
                    dist = harmonic_mean_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "cosine") {
                    dist = cosine_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "hassebrook") {
                    dist = hassebrook_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "jaccard") {
                    dist = jaccard_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "dice") {
                    dist = dice_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "fidelity") {
                    dist = fidelity_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "bhattacharyya") {
                    dist = bhattacharyya_internal(row_i.begin(), row_i.end(), row_j.begin(), unit, epsilon);
                } else if (method == "hellinger") {
                    dist = hellinger_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "matusita") {
                    dist = matusita_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "squared_chord") {
                    dist = squared_chord_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "squared_euclidean") {
                    dist = squared_euclidean_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "pearson") {
                    dist = pearson_internal(row_i.begin(), row_i.end(), row_j.begin(), epsilon);
                } else if (method == "neyman") {
                    dist = neyman_internal(row_i.begin(), row_i.end(), row_j.begin(), epsilon);
                } else if (method == "squared_chi") {
                    dist = squared_chi_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "prob_symm") {
                    dist = prob_symm_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "divergence") {
                    dist = divergence_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "clark") {
                    dist = clark_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "additive_symm") {
                    dist = additive_symm_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "kullback-leibler") {
                    dist = kullback_leibler_internal(row_i.begin(), row_i.end(), row_j.begin(), unit, epsilon);
                } else if (method == "jeffreys") {
                    dist = jeffreys_internal(row_i.begin(), row_i.end(), row_j.begin(), unit, epsilon);
                } else if (method == "k_divergence") {
                    dist = k_divergence_internal(row_i.begin(), row_i.end(), row_j.begin(), unit);
                } else if (method == "topsoe") {
                    dist = topsoe_internal(row_i.begin(), row_i.end(), row_j.begin(), unit);
                } else if (method == "jensen-shannon"){
                    dist = jensen_shannon_internal(row_i.begin(), row_i.end(), row_j.begin(), unit);
                } else if (method == "jensen_difference") {
                    dist = jensen_difference_internal(row_i.begin(), row_i.end(), row_j.begin(), unit);
                } else if (method == "taneja") {
                    dist = taneja_internal(row_i.begin(), row_i.end(), row_j.begin(), unit, epsilon);
                } else {
                    Rcpp::stop("Method not implemented for parallel execution.");
                }
                dist_matrix(i, j) = dist;
            }
        }
    }
};

Rcpp::NumericMatrix dist_many_many_parallel(Rcpp::NumericMatrix& dists1, Rcpp::NumericMatrix& dists2, Rcpp::String method, double p = NA_REAL, bool testNA = true, Rcpp::String unit = "log", double epsilon = 0.00001){
    int n1 = dists1.nrow();
    int n2 = dists2.nrow();
    Rcpp::NumericMatrix dist_matrix(n1, n2);
    std::string method_str(method);
    std::string unit_str(unit);

    ManyManyWorker worker(dists1, dists2, dist_matrix, method_str, p, unit_str, epsilon);
    RcppParallel::parallelFor(0, n1, worker);
    return dist_matrix;
}


// Worker and function for DistMatrixWithoutUnitMAT_parallel
struct DistMatrixWithoutUnitMATWorker : public RcppParallel::Worker {
    RcppParallel::RMatrix<double> dists;
    RcppParallel::RMatrix<double> dist_matrix;
    std::string method;
    double p;

    DistMatrixWithoutUnitMATWorker(Rcpp::NumericMatrix& dists,
                           Rcpp::NumericMatrix& dist_matrix,
                           std::string method,
                           double p)
        : dists(dists), dist_matrix(dist_matrix), method(method), p(p) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            for (std::size_t j = i; j < (std::size_t)dists.ncol(); ++j) {
                double dist = 0.0;
                auto col_i = dists.column(i);
                auto col_j = dists.column(j);

                if (method == "euclidean") {
                    dist = euclidean_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "minkowski") {
                    dist = minkowski_internal(col_i.begin(), col_i.end(), col_j.begin(), p);
                } else if (method == "manhattan") {
                    // manhattan is minkowski with p = 1
                    dist = minkowski_internal(col_i.begin(), col_i.end(), col_j.begin(), 1.0);
                } else if (method == "chebyshev") {
                    double max_dist = 0.0;
                    for(size_t k = 0; k < (size_t)col_i.length(); ++k) {
                        max_dist = std::max(max_dist, std::abs(col_i[k] - col_j[k]));
                    }
                    dist = max_dist;
                } else if (method == "lorentzian") {
                    dist = lorentzian_internal(col_i.begin(), col_i.end(), col_j.begin(), "log");
                } else if (method == "canberra") {
                    dist = canberra_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "sorensen") {
                    dist = sorensen_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "gower") {
                    dist = gower_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "soergel") {
                    dist = soergel_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "kulczynski_d") {
                    dist = kulczynski_d_internal(col_i.begin(), col_i.end(), col_j.begin(), 0.00001);
                } else if (method == "intersection") {
                    dist = intersection_dist_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "wavehedges") {
                    dist = wave_hedges_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "czekanowski") {
                    dist = czekanowski_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "motyka") {
                    dist = motyka_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "tanimoto") {
                    dist = tanimoto_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "ruzicka") {
                    dist = ruzicka_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "inner_product") {
                    dist = inner_product_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "harmonic_mean") {
                    dist = harmonic_mean_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "cosine") {
                    dist = cosine_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "hassebrook") {
                    dist = hassebrook_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "jaccard") {
                    dist = jaccard_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "dice") {
                    dist = dice_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "fidelity") {
                    dist = fidelity_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "hellinger") {
                    dist = hellinger_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "matusita") {
                    dist = matusita_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "squared_chord") {
                    dist = squared_chord_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "squared_euclidean") {
                    dist = squared_euclidean_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "pearson") {
                    dist = pearson_internal(col_i.begin(), col_i.end(), col_j.begin(), 0.00001);
                } else if (method == "neyman") {
                    dist = neyman_internal(col_i.begin(), col_i.end(), col_j.begin(), 0.00001);
                } else if (method == "squared_chi") {
                    dist = squared_chi_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "prob_symm") {
                    dist = prob_symm_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "divergence") {
                    dist = divergence_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "clark") {
                    dist = clark_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "additive_symm") {
                    dist = additive_symm_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else if (method == "kumar-johnson") {
                    dist = kumar_johnson_internal(col_i.begin(), col_i.end(), col_j.begin(), 0.00001);
                } else if (method == "avg") {
                    dist = avg_internal(col_i.begin(), col_i.end(), col_j.begin());
                } else {
                    Rcpp::stop("Method not implemented for parallel execution without unit.");
                }
                dist_matrix(i, j) = dist;
                dist_matrix(j, i) = dist;
            }
        }
    }
};

Rcpp::NumericMatrix DistMatrixWithoutUnitMAT_internal(Rcpp::NumericMatrix dists,
                                                    std::string method,
                                                    bool testNA,
                                                    Rcpp::Nullable<double> p = R_NilValue) {
    int n = dists.ncol();
    Rcpp::NumericMatrix dist_matrix(n, n);
    double p_val = 2.0;
    if (p.isNotNull()) {
        p_val = Rcpp::as<double>(p);
    }

    DistMatrixWithoutUnitMATWorker worker(dists, dist_matrix, method, p_val);
    RcppParallel::parallelFor(0, n, worker);

    return dist_matrix;
}

// Worker and function for DistMatrixWithUnitMAT_parallel
struct DistWorkerWithUnitMAT : public RcppParallel::Worker {
    RcppParallel::RMatrix<double> dists;
    RcppParallel::RMatrix<double> dist_matrix;
    std::string method;
    std::string unit;

    DistWorkerWithUnitMAT(Rcpp::NumericMatrix& dists,
                             Rcpp::NumericMatrix& dist_matrix,
                             std::string method,
                             std::string unit)
        : dists(dists), dist_matrix(dist_matrix), method(method), unit(unit) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            for (std::size_t j = i; j < (std::size_t)dists.ncol(); ++j) {
                double dist = 0.0;
                auto col_i = dists.column(i);
                auto col_j = dists.column(j);

                if (method == "lorentzian") {
                    dist = lorentzian_internal(col_i.begin(), col_i.end(), col_j.begin(), unit);
                } else if (method == "bhattacharyya") {
                    dist = bhattacharyya_internal(col_i.begin(), col_i.end(), col_j.begin(), unit, 0.00001);
                } else if (method == "kullback-leibler") {
                    dist = kullback_leibler_internal(col_i.begin(), col_i.end(), col_j.begin(), unit, 0.00001);
                } else if (method == "jeffreys") {
                    dist = jeffreys_internal(col_i.begin(), col_i.end(), col_j.begin(), unit, 0.00001);
                } else if (method == "k_divergence") {
                    dist = k_divergence_internal(col_i.begin(), col_i.end(), col_j.begin(), unit);
                } else if (method == "topsoe") {
                    dist = topsoe_internal(col_i.begin(), col_i.end(), col_j.begin(), unit);
                } else if (method == "jensen-shannon") {
                    dist = jensen_shannon_internal(col_i.begin(), col_i.end(), col_j.begin(), unit);
                } else if (method == "jensen_difference") {
                    dist = jensen_difference_internal(col_i.begin(), col_i.end(), col_j.begin(), unit);
                } else if (method == "taneja") {
                    dist = taneja_internal(col_i.begin(), col_i.end(), col_j.begin(), unit, 0.00001);
                } else {
                    Rcpp::stop("Method not implemented for parallel execution with unit.");
                }
                dist_matrix(i, j) = dist;
                dist_matrix(j, i) = dist;
            }
        }
    }
};

Rcpp::NumericMatrix DistMatrixWithUnitMAT_internal(Rcpp::NumericMatrix dists,
                                                    std::string method,
                                                    bool testNA,
                                                    std::string unit) {
    int n = dists.ncol();
    Rcpp::NumericMatrix dist_matrix(n, n);

    DistWorkerWithUnitMAT worker(dists, dist_matrix, method, unit);
    RcppParallel::parallelFor(0, n, worker);

    return dist_matrix;
}

// Worker and function for DistMatrixWithoutUnitDF_parallel
struct NoUnitDFDistWorker : public RcppParallel::Worker {
    RcppParallel::RMatrix<double> dists;
    RcppParallel::RMatrix<double> dist_matrix;
    std::string method;
    double p;

    NoUnitDFDistWorker(Rcpp::NumericMatrix& dists,
                       Rcpp::NumericMatrix& dists_matrix,
                       std::string method,
                       double p)
        : dists(dists), dist_matrix(dists_matrix), method(method), p(p) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            for (std::size_t j = i; j < (std::size_t)dists.nrow(); ++j) {
                double dist = 0.0;
                if (method == "euclidean") {
                    auto row_i = dists.row(i);
                    auto row_j = dists.row(j);
                    dist = euclidean_internal(row_i.begin(), row_i.end(), row_j.begin());
                } else if (method == "minkowski") {
                    auto row_i = dists.row(i);
                    auto row_j = dists.row(j);
                    dist = minkowski_internal(row_i.begin(), row_i.end(), row_j.begin(), p);
                } else if (method == "manhattan") {
                    auto row_i = dists.row(i);
                    auto row_j = dists.row(j);
                    // manhattan is minkowski with p = 1
                    dist = minkowski_internal(row_i.begin(), row_i.end(), row_j.begin(), 1.0);
                } else if (method == "chebyshev") {
                    auto row_i = dists.row(i);
                    auto row_j = dists.row(j);
                    double max_dist = 0.0;
                    for(size_t k = 0; k < (size_t)row_i.length(); ++k) {
                        max_dist = std::max(max_dist, std::abs(row_i[k] - row_j[k]));
                    }
                    dist = max_dist;
                } else if (method == "sorensen") {
                    dist = sorensen_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "gower") {
                    dist = gower_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "soergel") {
                    dist = soergel_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "kulczynski_d") {
                    dist = kulczynski_d_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin(), 0.00001);
                } else if (method == "canberra") {
                    dist = canberra_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "intersection") {
                    dist = intersection_dist_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "wavehedges") {
                    dist = wave_hedges_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "czekanowski") {
                    dist = czekanowski_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "motyka") {
                    dist = motyka_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "tanimoto") {
                    dist = tanimoto_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "ruzicka") {
                    dist = ruzicka_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "inner_product") {
                    dist = inner_product_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "harmonic_mean") {
                    dist = harmonic_mean_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "cosine") {
                    dist = cosine_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "hassebrook") {
                    dist = hassebrook_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "jaccard") {
                    dist = jaccard_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "dice") {
                    dist = dice_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "fidelity") {
                    dist = fidelity_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "hellinger") {
                    dist = hellinger_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "matusita") {
                    dist = matusita_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "squared_chord") {
                    dist = squared_chord_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "squared_euclidean") {
                    dist = squared_euclidean_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else if (method == "kumar-johnson") {
                    dist = kumar_johnson_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin(), 0.00001);
                } else if (method == "avg") {
                    dist = avg_internal(dists.row(i).begin(), dists.row(i).end(), dists.row(j).begin());
                } else {
                    Rcpp::stop("Method not implemented for parallel execution without unit.");
                }
                dist_matrix(i, j) = dist;
                dist_matrix(j, i) = dist;
            }
        }
    }
};


Rcpp::NumericMatrix DistMatrixWithoutUnitDF_internal(Rcpp::DataFrame distsDF,
                                                    std::string method,
                                                    bool testNA,
                                                    Rcpp::Nullable<double> p = R_NilValue) {
    Rcpp::NumericMatrix dists = as_matrix(distsDF);
    int n = dists.nrow();
    Rcpp::NumericMatrix dist_matrix(n, n);
    double p_val = 2.0;
    if (p.isNotNull()) {
        p_val = Rcpp::as<double>(p);
    }

    NoUnitDFDistWorker worker(dists, dist_matrix, method, p_val);
    RcppParallel::parallelFor(0, n, worker);

    return dist_matrix;
}

// Worker and function for DistMatrixWithUnitDF_parallel
struct UnitDFDistWorker : public RcppParallel::Worker {
    RcppParallel::RMatrix<double> dists;
    RcppParallel::RMatrix<double> dist_matrix;
    std::string method;
    std::string unit;

    UnitDFDistWorker(Rcpp::NumericMatrix& dists,
                     Rcpp::NumericMatrix& dist_matrix,
                     std::string method,
                     std::string unit)
        : dists(dists), dist_matrix(dist_matrix), method(method), unit(unit) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            for (std::size_t j = i; j < (std::size_t)dists.nrow(); ++j) {
                double dist = 0.0;
                auto row_i = dists.row(i);
                auto row_j = dists.row(j);
                if (method == "lorentzian") {
                    dist = lorentzian_internal(row_i.begin(), row_i.end(), row_j.begin(), unit);
                } else if (method == "bhattacharyya") {
                    dist = bhattacharyya_internal(row_i.begin(), row_i.end(), row_j.begin(), unit, 0.00001);
                } else if (method == "kullback-leibler") {
                    dist = kullback_leibler_internal(row_i.begin(), row_i.end(), row_j.begin(), unit, 0.00001);
                } else if (method == "jeffreys") {
                    dist = jeffreys_internal(row_i.begin(), row_i.end(), row_j.begin(), unit, 0.00001);
                } else if (method == "k_divergence") {
                    dist = k_divergence_internal(row_i.begin(), row_i.end(), row_j.begin(), unit);
                } else if (method == "topsoe") {
                    dist = topsoe_internal(row_i.begin(), row_i.end(), row_j.begin(), unit);
                } else if (method == "jensen-shannon") {
                    dist = jensen_shannon_internal(row_i.begin(), row_i.end(), row_j.begin(), unit);
                } else if (method == "jensen_difference") {
                    dist = jensen_difference_internal(row_i.begin(), row_i.end(), row_j.begin(), unit);
                } else if (method == "taneja") {
                    dist = taneja_internal(row_i.begin(), row_i.end(), row_j.begin(), unit, 0.00001);
                } else {
                    Rcpp::stop("Method not implemented for parallel execution with unit.");
                }
                dist_matrix(i, j) = dist;
                dist_matrix(j, i) = dist;
            }
        }
    }
};

Rcpp::NumericMatrix DistMatrixWithUnitDF_internal(Rcpp::DataFrame distsDF,
                                                    std::string method,
                                                    bool testNA,
                                                    std::string unit) {
    Rcpp::NumericMatrix dists = as_matrix(distsDF);
    int n = dists.nrow();
    Rcpp::NumericMatrix dist_matrix(n, n);

    UnitDFDistWorker worker(dists, dist_matrix, method, unit);
    RcppParallel::parallelFor(0, n, worker);

    return dist_matrix;
}