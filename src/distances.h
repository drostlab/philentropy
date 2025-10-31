//  Part of the philentropy package
//
//  Copyright (C) 2015-2020 Hajk-Georg Drost
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  A copy of the GNU General Public License is available at
//  http://www.r-project.org/Licenses/
// 

#ifndef philentropy_Distances_H
#define philentropy_Distances_H philentropy_Distances_H

// http://stackoverflow.com/questions/23527719/calling-a-rcpp-function-from-another-rcpp-function-while-building-an-r-package


// [[Rcpp::plugins(cpp11)]]


#include <Rcpp.h> 
#include <math.h>
#include <iostream>
#include "utils.h"

// [[Rcpp::export]]
double custom_log2(const double& x ){
        if (x == 0.0){
          return NAN;
        } else {
          return log(x)/log(2.0);
        }
}

// [[Rcpp::export]]
double custom_log10(const double& x ){
  if (x == 0.0){
    return NAN;
  } else {
    return log(x)/log(10.0);
  }
}

//' @title Euclidean distance (lowlevel function)
//' @description The lowlevel function for computing the euclidean distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' euclidean(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double euclidean(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return euclidean_internal(P.begin(), P.end(), Q.begin());
}

//' @title Manhattan distance (lowlevel function)
//' @description The lowlevel function for computing the manhattan distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' manhattan(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double manhattan(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return minkowski_internal(P.begin(), P.end(), Q.begin(), 1.0);
}



//' @title Minkowski distance (lowlevel function)
//' @description The lowlevel function for computing the minkowski distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param n index for the minkowski exponent.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' minkowski(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), n = 2, testNA = FALSE)
//' @export
// [[Rcpp::export]]
double minkowski(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, double n, bool testNA = true){
        if (testNA) check_na(P, Q);
        validate_p_parameter("minkowski", n);
        return minkowski_internal(P.begin(), P.end(), Q.begin(), n);
}


//' @title Chebyshev distance (lowlevel function)
//' @description The lowlevel function for computing the chebyshev distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' chebyshev(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double chebyshev(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return chebyshev_internal(P.begin(), P.end(), Q.begin());
}


//' @title Sorensen distance (lowlevel function)
//' @description The lowlevel function for computing the sorensen distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' sorensen(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double sorensen(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return sorensen_internal(P.begin(), P.end(), Q.begin());
}



//' @title Gower distance (lowlevel function)
//' @description The lowlevel function for computing the gower distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' gower(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double gower(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return gower_internal(P.begin(), P.end(), Q.begin());
}



//' @title Soergel distance (lowlevel function)
//' @description The lowlevel function for computing the soergel distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' soergel(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double soergel(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return soergel_internal(P.begin(), P.end(), Q.begin());
}


//' @title Kulczynski_d distance (lowlevel function)
//' @description The lowlevel function for computing the kulczynski_d distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
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
//' @author Hajk-Georg Drost
//' @examples
//' kulczynski_d(P = 1:10/sum(1:10), Q = 20:29/sum(20:29),
//'     testNA = FALSE, epsilon = 0.00001)
//' @export
// [[Rcpp::export]]
double kulczynski_d(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true, double epsilon = 0.00001){
        if (testNA) check_na(P, Q);
        return kulczynski_d_internal(P.begin(), P.end(), Q.begin(), epsilon);
}


//' @title Canberra distance (lowlevel function)
//' @description The lowlevel function for computing the canberra distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' canberra(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double canberra(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return canberra_internal(P.begin(), P.end(), Q.begin());
}


//' @title Lorentzian distance (lowlevel function)
//' @description The low-level function for computing the lorentzian distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @param unit type of \code{log} function. Option are 
//' \itemize{
//' \item \code{unit = "log"}
//' \item \code{unit = "log2"}
//' \item \code{unit = "log10"}   
//' }
//' @author Hajk-Georg Drost
//' @examples
//' lorentzian(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE, unit = "log2")
//' @export
// [[Rcpp::export]]
double lorentzian(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true, const Rcpp::String unit = "log"){
        if (testNA) check_na(P, Q);
        std::string unit_str(unit);
        return lorentzian_internal(P.begin(), P.end(), Q.begin(), unit_str);
}


//' @title Intersection distance (lowlevel function)
//' @description The lowlevel function for computing the intersection_dist distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' intersection_dist(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double intersection_dist(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return intersection_dist_internal(P.begin(), P.end(), Q.begin());
}


//' @title Wave hedges distance (lowlevel function)
//' @description The lowlevel function for computing the wave_hedges distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' wave_hedges(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double wave_hedges(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return wave_hedges_internal(P.begin(), P.end(), Q.begin());
}


//' @title Czekanowski distance (lowlevel function)
//' @description The lowlevel function for computing the czekanowski distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' czekanowski(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double czekanowski(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return czekanowski_internal(P.begin(), P.end(), Q.begin());
}



//' @title Motyka distance (lowlevel function)
//' @description The lowlevel function for computing the motyka distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' motyka(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double motyka(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return motyka_internal(P.begin(), P.end(), Q.begin());
}


//' @title Tanimoto distance (lowlevel function)
//' @description The lowlevel function for computing the tanimoto distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' tanimoto(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double tanimoto(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return tanimoto_internal(P.begin(), P.end(), Q.begin());
}


//' @title Ruzicka distance (lowlevel function)
//' @description The lowlevel function for computing the ruzicka distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' ruzicka(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double ruzicka(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return ruzicka_internal(P.begin(), P.end(), Q.begin());
}


//' @title Inner product distance (lowlevel function)
//' @description The lowlevel function for computing the inner_product distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' inner_product(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double inner_product(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return inner_product_internal(P.begin(), P.end(), Q.begin());
}


//' @title Harmonic mean distance (lowlevel function)
//' @description The lowlevel function for computing the harmonic_mean_dist distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' harmonic_mean_dist(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double harmonic_mean_dist(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return harmonic_mean_internal(P.begin(), P.end(), Q.begin());
}


//' @title Cosine distance (lowlevel function)
//' @description The lowlevel function for computing the cosine_dist distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' cosine_dist(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double cosine_dist(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return cosine_internal(P.begin(), P.end(), Q.begin());
}


//' @title Kumar hassebrook distance (lowlevel function)
//' @description The lowlevel function for computing the kumar_hassebrook distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' kumar_hassebrook(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double kumar_hassebrook(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return hassebrook_internal(P.begin(), P.end(), Q.begin());
}


//' @title Jaccard distance (lowlevel function)
//' @description The lowlevel function for computing the jaccard distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' jaccard(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double jaccard(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return jaccard_internal(P.begin(), P.end(), Q.begin());
}


//' @title Dice distance (lowlevel function)
//' @description The lowlevel function for computing the dice_dist distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' dice_dist(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double dice_dist(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return dice_internal(P.begin(), P.end(), Q.begin());
}


//' @title Fidelity distance (lowlevel function)
//' @description The lowlevel function for computing the fidelity distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' fidelity(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double fidelity(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return fidelity_internal(P.begin(), P.end(), Q.begin());
}


//' @title Bhattacharyya distance (lowlevel function)
//' @description The lowlevel function for computing the bhattacharyya distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @param unit type of \code{log} function. Option are 
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
//' @author Hajk-Georg Drost
//' @examples
//' bhattacharyya(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE,
//'  unit = "log2", epsilon = 0.00001)
//' @export
// [[Rcpp::export]]
double bhattacharyya(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true, const Rcpp::String unit = "log", double epsilon = 0.00001){
        if (testNA) check_na(P, Q);
        std::string unit_str(unit);
        return bhattacharyya_internal(P.begin(), P.end(), Q.begin(), unit_str, epsilon);
}


//' @title Hellinger distance (lowlevel function)
//' @description The lowlevel function for computing the hellinger distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' hellinger(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double hellinger(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return hellinger_internal(P.begin(), P.end(), Q.begin());
}

//' @title Matusita distance (lowlevel function)
//' @description The lowlevel function for computing the matusita distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' matusita(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double matusita(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return matusita_internal(P.begin(), P.end(), Q.begin());
}

//' @title Squared chord distance (lowlevel function)
//' @description The lowlevel function for computing the squared_chord distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' squared_chord(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double squared_chord(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return squared_chord_internal(P.begin(), P.end(), Q.begin());
}



//' @title Squared euclidean distance (lowlevel function)
//' @description The lowlevel function for computing the squared_euclidean distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' squared_euclidean(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double squared_euclidean(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return squared_euclidean_internal(P.begin(), P.end(), Q.begin());
}


//' @title Pearson chi-squared distance (lowlevel function)
//' @description The lowlevel function for computing the pearson_chi_sq distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
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
//' @author Hajk-Georg Drost
//' @examples
//' pearson_chi_sq(P = 1:10/sum(1:10), Q = 20:29/sum(20:29),
//'  testNA = FALSE, epsilon = 0.00001)
//' @export
// [[Rcpp::export]]
double pearson_chi_sq(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true, double epsilon = 0.00001){
        if (testNA) check_na(P, Q);
        return pearson_internal(P.begin(), P.end(), Q.begin(), epsilon);
}


//' @title Neyman chi-squared distance (lowlevel function)
//' @description The lowlevel function for computing the neyman_chi_sq distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
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
//' @author Hajk-Georg Drost
//' @examples
//' neyman_chi_sq(P = 1:10/sum(1:10), Q = 20:29/sum(20:29),
//'  testNA = FALSE, epsilon = 0.00001)
//' @export
// [[Rcpp::export]]
double neyman_chi_sq(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true, double epsilon = 0.00001){
        if (testNA) check_na(P, Q);
        return neyman_internal(P.begin(), P.end(), Q.begin(), epsilon);
}


//' @title Squared chi-squared distance (lowlevel function)
//' @description The lowlevel function for computing the squared_chi_sq distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' squared_chi_sq(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double squared_chi_sq(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return squared_chi_internal(P.begin(), P.end(), Q.begin());
}


//' @title Probability symmetric chi-squared distance (lowlevel function)
//' @description The lowlevel function for computing the prob_symm_chi_sq distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' prob_symm_chi_sq(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double prob_symm_chi_sq(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return prob_symm_internal(P.begin(), P.end(), Q.begin());
}




//' @title Divergence squared distance (lowlevel function)
//' @description The lowlevel function for computing the divergence_sq distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' divergence_sq(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double divergence_sq(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return divergence_internal(P.begin(), P.end(), Q.begin());
}



//' @title Clark squared distance (lowlevel function)
//' @description The lowlevel function for computing the clark_sq distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' clark_sq(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double clark_sq(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return clark_internal(P.begin(), P.end(), Q.begin());
}


//' @title Additive symmetric chi-squared distance (lowlevel function)
//' @description The lowlevel function for computing the additive_symm_chi_sq distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' additive_symm_chi_sq(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double additive_symm_chi_sq(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return additive_symm_internal(P.begin(), P.end(), Q.begin());
}



//' @title kullback-Leibler distance (lowlevel function)
//' @description The lowlevel function for computing the kullback_leibler_distance distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
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
//' @author Hajk-Georg Drost
//' @examples
//' kullback_leibler_distance(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE,
//'  unit = "log2", epsilon = 0.00001)
//' @export
// [[Rcpp::export]]
double kullback_leibler_distance(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true, const Rcpp::String unit = "log", double epsilon = 0.00001){
        if (testNA) check_na(P, Q);
        std::string unit_str(unit);
        return kullback_leibler_internal(P.begin(), P.end(), Q.begin(), unit_str, epsilon);
}

//' @title Jeffreys distance (lowlevel function)
//' @description The lowlevel function for computing the jeffreys distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
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
//' @author Hajk-Georg Drost
//' @examples
//' jeffreys(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE,
//'  unit = "log2", epsilon = 0.00001)
//' @export
// [[Rcpp::export]]
double jeffreys(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true, const Rcpp::String unit = "log", double epsilon = 0.00001){
        if (testNA) check_na(P, Q);
        std::string unit_str(unit);
        return jeffreys_internal(P.begin(), P.end(), Q.begin(), unit_str, epsilon);
}


//' @title K-Divergence (lowlevel function)
//' @description The lowlevel function for computing the k_divergence distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @param unit type of \code{log} function. Option are 
//' \itemize{
//' \item \code{unit = "log"}
//' \item \code{unit = "log2"}
//' \item \code{unit = "log10"}   
//' }
//' @author Hajk-Georg Drost
//' @examples
//' k_divergence(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE, unit = "log2")
//' @export
// [[Rcpp::export]]
double k_divergence(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true, const Rcpp::String unit = "log"){
        if (testNA) check_na(P, Q);
        std::string unit_str(unit);
        return k_divergence_internal(P.begin(), P.end(), Q.begin(), unit_str);
}


//' @title Topsoe distance (lowlevel function)
//' @description The lowlevel function for computing the topsoe distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @param unit type of \code{log} function. Option are 
//' \itemize{
//' \item \code{unit = "log"}
//' \item \code{unit = "log2"}
//' \item \code{unit = "log10"}   
//' }
//' @author Hajk-Georg Drost
//' @examples
//' topsoe(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE, unit = "log2")
//' @export
// [[Rcpp::export]]
double topsoe(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true, const Rcpp::String unit = "log"){
        if (testNA) check_na(P, Q);
        std::string unit_str(unit);
        return topsoe_internal(P.begin(), P.end(), Q.begin(), unit_str);
}


//' @title Jensen-Shannon distance (lowlevel function)
//' @description The lowlevel function for computing the jensen_shannon distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @param unit type of \code{log} function. Option are 
//' \itemize{
//' \item \code{unit = "log"}
//' \item \code{unit = "log2"}
//' \item \code{unit = "log10"}   
//' }
//' @author Hajk-Georg Drost
//' @examples
//' jensen_shannon(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE, unit = "log2")
//' @export
// [[Rcpp::export]]
double jensen_shannon(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true, const Rcpp::String unit = "log"){
        if (testNA) check_na(P, Q);
        std::string unit_str(unit);
        return jensen_shannon_internal(P.begin(), P.end(), Q.begin(), unit_str);
}


//' @title Jensen difference (lowlevel function)
//' @description The lowlevel function for computing the jensen_difference distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @param unit type of \code{log} function. Option are 
//' \itemize{
//' \item \code{unit = "log"}
//' \item \code{unit = "log2"}
//' \item \code{unit = "log10"}   
//' }
//' @author Hajk-Georg Drost
//' @examples
//' jensen_difference(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE, unit = "log2")
//' @export
// [[Rcpp::export]]
double jensen_difference(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true, const Rcpp::String unit = "log"){
        if (testNA) check_na(P, Q);
        std::string unit_str(unit);
        return jensen_difference_internal(P.begin(), P.end(), Q.begin(), unit_str);
}




//' @title Taneja difference (lowlevel function)
//' @description The lowlevel function for computing the taneja distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
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
//' @author Hajk-Georg Drost
//' @examples
//' taneja(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE,
//'  unit = "log2", epsilon = 0.00001)
//' @export
// [[Rcpp::export]]
double taneja(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true, const Rcpp::String unit = "log", double epsilon = 0.00001){
        if (testNA) check_na(P, Q);
        std::string unit_str(unit);
        return taneja_internal(P.begin(), P.end(), Q.begin(), unit_str, epsilon);
}


//' @title Kumar-Johnson distance (lowlevel function)
//' @description The lowlevel function for computing the kumar_johnson distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
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
//' @author Hajk-Georg Drost
//' @examples
//' kumar_johnson(P = 1:10/sum(1:10), Q = 20:29/sum(20:29),
//'  testNA = FALSE, epsilon = 0.00001)
//' @export
// [[Rcpp::export]]
double kumar_johnson(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true, double epsilon = 0.00001){
        if (testNA) check_na(P, Q);
        return kumar_johnson_internal(P.begin(), P.end(), Q.begin(), epsilon);
}



//' @title AVG distance (lowlevel function)
//' @description The lowlevel function for computing the avg distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' avg(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double avg(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA = true){
        if (testNA) check_na(P, Q);
        return avg_internal(P.begin(), P.end(), Q.begin());
}

#endif // philentropy_Distances_H
