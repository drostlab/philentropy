

#ifndef philentropy_UTILS_H
#define philentropy_UTILS_H philentropy_UTILS_H

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <set>
#include <string>

Rcpp::NumericMatrix as_matrix(Rcpp::DataFrame x);

Rcpp::DataFrame as_data_frame(Rcpp::NumericMatrix mat);

SEXP sum_rcpp(SEXP vec);

SEXP est_prob_empirical(SEXP CountVec);

inline void check_na(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q) {
    if (P.size() != Q.size())
        Rcpp::stop("Input vectors do not have the same length.");

    for (R_xlen_t i = 0; i < P.size(); ++i) {
        if (Rcpp::NumericVector::is_na(P[i]) || Rcpp::NumericVector::is_na(Q[i]))
            Rcpp::stop("Input vectors contain NA values.");
    }
}


// validate 'p' for methods that require it before entering parallel regions
inline void validate_p_parameter(const std::string& method, double p) {
    const std::set<std::string> p_methods = {"minkowski"};
    if (p_methods.count(method) && std::isnan(p)) {
        Rcpp::stop("Please specify the 'p' parameter for the '" + method + "' distance.");
    }
}

// C++ helper to resolve number of threads
inline int get_num_threads_cpp(Rcpp::Nullable<int> num_threads) {
    if (num_threads.isNotNull()) {
        return Rcpp::as<int>(num_threads);
    } else {
        char* env_var = std::getenv("RCPP_PARALLEL_NUM_THREADS");
        if (env_var) {
            return std::atoi(env_var);
        } else {
            return 2;
        }
    }
}

#endif // philentropy_UTILS_H
