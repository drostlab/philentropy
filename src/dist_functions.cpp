// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::interfaces(r, cpp)]] 

#include <Rcpp.h>
#include "utils.h"
#include "dist_dispatch.h"
#include "InformationTheory.h"
#include "correlation.h"
#include <RcppParallel.h>

//' @title Distances and Similarities between Two Probability Density Functions
//' @name dist_one_one
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
        std::string unit_str(unit);
        std::string method_str(method);
        if (method_str == "minkowski") {
            if (p.isNotNull()) {
                return minkowski_internal(P.begin(), P.end(), Q.begin(), Rcpp::as<double>(p));
            } else {
                Rcpp::stop("Please specify p for the Minkowski distance.");
            }
        }
        return dispatch_dist_internal(P.begin(), P.end(), Q.begin(), method_str, unit_str, epsilon, NAN);
}

struct OneManyWorker : public RcppParallel::Worker {
    const RcppParallel::RVector<double> P;
    RcppParallel::RMatrix<double> dists;
    RcppParallel::RVector<double> dist_values;

    const std::string method;
    const bool testNA;
    const std::string unit;
    const double epsilon;
    const double p;

    OneManyWorker(const Rcpp::NumericVector& P_in,
                  const Rcpp::NumericMatrix& dists_in,
                  Rcpp::NumericVector& dist_values_in,
                  const Rcpp::String& method_in,
                  const bool& testNA_in,
                  const Rcpp::String& unit_in,
                  const double& epsilon_in,
                  const double p_in)
        : P(P_in), dists(dists_in), dist_values(dist_values_in),
          method(method_in), testNA(testNA_in),
          unit(unit_in), epsilon(epsilon_in), p(p_in) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            auto Q = dists.row(i);
            dist_values[i] = dispatch_dist_internal(P.begin(), P.end(), Q.begin(), method, unit, epsilon, p);
        }
    }
};

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
//' @param num_threads an integer specifying the number of threads to be used for parallel computations. Default is taken from the `RCPP_PARALLEL_NUM_THREADS` environment variable, or `2` if not set.
//' @return A vector of distance values
//' @examples
//'   set.seed(2020-08-20)
//'   P <- 1:10 / sum(1:10)
//'   M <- t(replicate(100, sample(1:10, size = 10) / 55))
//'   dist_one_many(P, M, method = "euclidean", testNA = FALSE)
//' @export
// [[Rcpp::export(name = "dist_one_many")]]
Rcpp::NumericVector dist_one_many_cpp(const Rcpp::NumericVector& P, Rcpp::NumericMatrix dists, Rcpp::String method, Rcpp::Nullable<double> p = R_NilValue, bool testNA = true, Rcpp::String unit = "log", double epsilon = 0.00001, Rcpp::Nullable<int> num_threads = R_NilValue) {
    std::string method_str(method);
    int nrows = dists.nrow();
    Rcpp::NumericVector dist_values(nrows);
    int n_threads = get_num_threads_cpp(num_threads);
    double p_val = p.isNotNull() ? Rcpp::as<double>(p) : NAN;

    validate_p_parameter(method_str, p_val);
    OneManyWorker worker(P, dists, dist_values, method, testNA, unit, epsilon, p_val);
    RcppParallel::parallelFor(0, nrows, worker, 1, n_threads);

    return dist_values;
}

struct ManyManyWorker : public RcppParallel::Worker {
    RcppParallel::RMatrix<double> dists1;
    RcppParallel::RMatrix<double> dists2;
    RcppParallel::RMatrix<double> dist_matrix;
    std::string method;
    std::string unit;
    double epsilon;
    double p;

    ManyManyWorker(Rcpp::NumericMatrix& dists1,
                   Rcpp::NumericMatrix& dists2,
                   Rcpp::NumericMatrix& dist_matrix,
                   std::string method,
                   std::string unit,
                   double epsilon,
                   double p)
        : dists1(dists1), dists2(dists2), dist_matrix(dist_matrix),
          method(method), unit(unit), epsilon(epsilon), p(p) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            for (std::size_t j = 0; j < (std::size_t)dists2.nrow(); ++j) {
                auto row_i = dists1.row(i);
                auto row_j = dists2.row(j);
                double dist = dispatch_dist_internal(row_i.begin(), row_i.end(), row_j.begin(),
                                                     method, unit, epsilon, p);
                dist_matrix(i, j) = dist;
            }
        }
    }
};

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
//' @param num_threads an integer specifying the number of threads to be used for parallel computations. Default is taken from the `RCPP_PARALLEL_NUM_THREADS` environment variable, or `2` if not set.
//' @return A matrix of distance values
//' @examples
//'   set.seed(2020-08-20)
//'   M1 <- t(replicate(10, sample(1:10, size = 10) / 55))
//'   M2 <- t(replicate(10, sample(1:10, size = 10) / 55))
//'   result <- dist_many_many(M1, M2, method = "euclidean", testNA = FALSE)
//' @export
// [[Rcpp::export(name = "dist_many_many")]]
Rcpp::NumericMatrix dist_many_many_cpp(Rcpp::NumericMatrix& dists1, Rcpp::NumericMatrix& dists2, Rcpp::String method, Rcpp::Nullable<double> p = R_NilValue, bool testNA = true, Rcpp::String unit = "log", double epsilon = 0.00001, Rcpp::Nullable<int> num_threads = R_NilValue) {
    int n1 = dists1.nrow();
    int n2 = dists2.nrow();
    Rcpp::NumericMatrix dist_matrix(n1, n2);
    std::string method_str(method);
    std::string unit_str(unit);
    int n_threads = get_num_threads_cpp(num_threads);
    double p_val = p.isNotNull() ? Rcpp::as<double>(p) : NAN;

    validate_p_parameter(method_str, p_val);
    ManyManyWorker worker(dists1, dists2, dist_matrix, method_str, unit_str, epsilon, p_val);
    RcppParallel::parallelFor(0, n1, worker, 1, n_threads);
    return dist_matrix;
}
