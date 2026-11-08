// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include "utils.h"
#include "dist_dispatch.h"
#include <RcppParallel.h>

struct DistMatrixWorker : public RcppParallel::Worker {
    RcppParallel::RMatrix<double> dists;
    RcppParallel::RMatrix<double> dist_matrix;
    std::string method;
    double epsilon;
    double p;

    DistMatrixWorker(Rcpp::NumericMatrix& dists,
                            Rcpp::NumericMatrix& dist_matrix,
                            std::string method,
                            double epsilon,
                            double p)
        : dists(dists), dist_matrix(dist_matrix), method(method), epsilon(epsilon), p(p) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            for (std::size_t j = i; j < (std::size_t)dists.ncol(); ++j) {
                double dist = 0.0;
                auto col_i = dists.column(i);
                auto col_j = dists.column(j);
                dist = dispatch_dist_internal(col_i.begin(), col_i.end(), col_j.begin(),
                                              method, "log", epsilon, p);
                dist_matrix(i, j) = dist;
                dist_matrix(j, i) = dist;
            }
        }
    }
};

Rcpp::NumericMatrix DistMatrixWithoutUnitMAT_internal(Rcpp::NumericMatrix dists,
                                                      std::string method,
                                                      bool testNA,
                                                      double epsilon,
                                                      Rcpp::Nullable<double> p,
                                                      Rcpp::Nullable<int> num_threads) {
    int n = dists.ncol();
    Rcpp::NumericMatrix dist_matrix(n, n);
    int n_threads = get_num_threads_cpp(num_threads);

    double p_val = NAN;
    if (p.isNotNull()) p_val = Rcpp::as<double>(p);
    validate_p_parameter(method, p_val);
    DistMatrixWorker worker(dists, dist_matrix, method, epsilon, p_val);
    RcppParallel::parallelFor(0, n, worker, 1, n_threads);

    return dist_matrix;
}

struct DistMatrixWorkerWithUnit : public RcppParallel::Worker {
    RcppParallel::RMatrix<double> dists;
    RcppParallel::RMatrix<double> dist_matrix;
    std::string method;
    std::string unit;
    double epsilon;

    DistMatrixWorkerWithUnit(Rcpp::NumericMatrix& dists,
                             Rcpp::NumericMatrix& dist_matrix,
                             std::string method,
                             std::string unit,
                             double epsilon)
        : dists(dists), dist_matrix(dist_matrix), method(method), unit(unit), epsilon(epsilon) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            for (std::size_t j = i; j < (std::size_t)dists.ncol(); ++j) {
                double dist = 0.0;
                auto col_i = dists.column(i);
                auto col_j = dists.column(j);
                dist = dispatch_dist_internal(col_i.begin(), col_i.end(), col_j.begin(), method, unit, epsilon, NAN);
                dist_matrix(i, j) = dist;
                dist_matrix(j, i) = dist;
            }
        }
    }
};

Rcpp::NumericMatrix DistMatrixWithUnitMAT_internal(Rcpp::NumericMatrix dists,
                                                   std::string method,
                                                   bool testNA,
                                                   double epsilon,
                                                   std::string unit,
                                                   Rcpp::Nullable<int> num_threads) {
    int n = dists.ncol();
    Rcpp::NumericMatrix dist_matrix(n, n);
    int n_threads = get_num_threads_cpp(num_threads);

    DistMatrixWorkerWithUnit worker(dists, dist_matrix, method, unit, epsilon);
    RcppParallel::parallelFor(0, n, worker, 1, n_threads);

    return dist_matrix;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix distance_cpp(Rcpp::NumericMatrix x,
                                 std::string method,
                                 Rcpp::Nullable<double> p,
                                 bool test_na,
                                 std::string unit,
                                 double epsilon,
                                 Rcpp::Nullable<int> num_threads) {

    // Define groups of methods
    const std::set<std::string> unit_methods = {
      "lorentzian", "bhattacharyya", "kullback-leibler", "jeffreys",
      "k_divergence", "topsoe", "jensen-shannon", "jensen_difference", "taneja"
    };

    if (unit_methods.count(method)) {
        return DistMatrixWithUnitMAT_internal(x, method, test_na, epsilon, unit, num_threads);
    } else if (method == "minkowski") {
        if (!p.isNotNull()) {
            Rcpp::stop("Please specify p for the Minkowski distance.");
        }
        return DistMatrixWithoutUnitMAT_internal(x, method, test_na, epsilon, p, num_threads);
    } else if (method == "non-intersection") {
        Rcpp::NumericMatrix intersection_matrix = DistMatrixWithoutUnitMAT_internal(x, "intersection", test_na, epsilon, p, num_threads);
        return 1.0 - intersection_matrix;
    } else if (method == "kulczynski_s") {
        Rcpp::NumericMatrix kulczynski_d_matrix = DistMatrixWithoutUnitMAT_internal(x, "kulczynski_d", test_na, epsilon, p, num_threads);
        // Element-wise division, handling potential division by zero
        for (int i = 0; i < kulczynski_d_matrix.nrow(); ++i) {
            for (int j = 0; j < kulczynski_d_matrix.ncol(); ++j) {
                if (kulczynski_d_matrix(i, j) != 0) {
                    kulczynski_d_matrix(i, j) = 1.0 / kulczynski_d_matrix(i, j);
                } else {
                    kulczynski_d_matrix(i, j) = R_PosInf; 
                }
            }
        }
        return kulczynski_d_matrix;
    } else {
        return DistMatrixWithoutUnitMAT_internal(x, method, test_na, epsilon, p, num_threads);
    }
}
