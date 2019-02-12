//  Part of the philentropy package
//
//  Copyright (C) 2015-2018 Hajk-Georg Drost
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
double euclidean(const Rcpp::NumericVector& P,const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0.0;
        double diff  = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }

        if (testNA){
                for (int i = 0; i < P_len; i++){     
                        if ((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        diff = fabs(P[i] - Q[i]);
                        dist += diff * diff;
                }
        } else {
              for (int i = 0; i < P_len; i++){
                      diff = fabs(P[i] - Q[i]);
                      dist += diff * diff;
               }  
        }
        return sqrt(dist);
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
double manhattan(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0.0;
        double diff  = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if (testNA){
                for (int i = 0; i < P_len; i++){
                        if ((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        diff = fabs(P[i] - Q[i]);
                        dist += diff;
                }   
        } else {
                for (int i = 0; i < P_len; i++){
                        diff = fabs(P[i] - Q[i]);
                        dist += diff;
                }  
        }
        return dist;
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
double minkowski(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, double n, bool testNA){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0.0;
        double diff  = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if (testNA){
                for (int i = 0; i < P_len; i++){
                        if ((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        diff = pow(fabs(P[i] - Q[i]), n);
                        dist += diff;
                }          
        } else {
                for (int i = 0; i < P_len; i++){
                        diff = pow(fabs(P[i] - Q[i]), n);
                        dist += diff;
                }
        }
        return pow(dist ,  1.0/n);
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
double chebyshev(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0.0;
        double diff  = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if (testNA){
                for (int i = 0; i < P_len; i++){
                        if ((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        diff = fabs(P[i] - Q[i]);
                        if (diff > dist)
                                dist = diff;
                }      
        } else {
               for (int i = 0; i < P_len; i++){
                        diff = fabs(P[i] - Q[i]);
                        if (diff > dist)
                                dist = diff;
                }   
        }
        return dist;
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
double sorensen(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist1 = 0.0;
        double dist2 = 0.0;
        double diff  = 0.0;
        double sum   = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if (testNA){
                for (int i = 0; i < P_len; i++){
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        diff = fabs(P[i] - Q[i]);
                        sum  = P[i] + Q[i];
                        dist1 += diff;
                        dist2 += sum;
                }  
        } else {
                for (int i = 0; i < P_len; i++){
                        diff = fabs(P[i] - Q[i]);
                        sum  = P[i] + Q[i];
                        dist1 += diff;
                        dist2 += sum;
                }
        }
        
        if (dist2 == 0.0) {
                return NAN;
        } else {
                return dist1/dist2;
        }
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
double gower(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double diff  = 0.0;
        double dist  = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if (P_len == 0){
                Rcpp::stop("One of the input vectors has length 0 and cannot be processed!");
        }
        
        if(testNA){
                for (int i = 0; i < P_len; i++){
                        if ((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        diff = fabs(P[i] - Q[i]);
                        dist += diff;
                }   
        } else {
                for (int i = 0; i < P_len; i++){
                        diff = fabs(P[i] - Q[i]);
                        dist += diff;
                } 
        }
        return ((1.0/P_len) * dist);    
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
double soergel(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double diff       = 0.0;
        double dist1      = 0.0;
        double dist2      = 0.0;
        double max_point  = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if (testNA){
                for (int i = 0; i < P_len; i++){
                        if ((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        diff      = fabs(P[i] - Q[i]);
                        if (P[i] >= Q[i]){
                                max_point = P[i];
                        } else {
                                max_point = Q[i];
                        }
                        dist1 += diff;
                        dist2 += max_point;
                }   
        } else {
                for (int i = 0; i < P_len; i++){
                        diff      = fabs(P[i] - Q[i]);
                        if (P[i] >= Q[i]){
                                max_point = P[i];
                        } else {
                                max_point = Q[i];
                        }
                        dist1 += diff;
                        dist2 += max_point;
                } 
        }
        
        if (dist2 == 0.0) {
                return 0;
        } else {
                return dist1/dist2;         
        }
}


//' @title Kulczynski_d distance (lowlevel function)
//' @description The lowlevel function for computing the kulczynski_d distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' kulczynski_d(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double kulczynski_d(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double diff       = 0.0;
        double dist1      = 0.0;
        double dist2      = 0.0;
        double min_point  = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if (testNA){
                for (int i = 0; i < P_len; i++){
                        if ((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        diff      = fabs(P[i] - Q[i]);
                        if (P[i] <= Q[i]){
                                min_point = P[i];
                        } else {
                                min_point = Q[i];
                        }
                        dist1 += diff;
                        if (min_point == 0.0){
                                dist2 += 0.00001;
                        } else {
                                dist2 += min_point;
                        }     
                }
        } else {
                for (int i = 0; i < P_len; i++){
                        diff      = fabs(P[i] - Q[i]);
                        if (P[i] <= Q[i]){
                                min_point = P[i];
                        } else {
                                min_point = Q[i];
                        }
                        dist1 += diff;
                        if (min_point == 0.0){
                                dist2 += 0.00001;
                        } else {
                                dist2 += min_point;
                        }     
                }
        }
        
        if (dist2 == 0.0) {
                return NAN;
        } else {
                return dist1/dist2;         
        }     
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
double canberra(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0.0;
        double diff  = 0.0;
        double sum   = 0.0;
        
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if (testNA){
                for (int i = 0; i < P_len; i++){
                        if ((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        diff = fabs(P[i] - Q[i]);
                        sum  = P[i] + Q[i];
                        // replace 0/0 or x/0 or 0/x by 0 according to Sung-Hyuk Cha (2007)
                        if ((diff == 0.0) || (sum == 0.0)){
                                dist += 0.0;
                        } else {
                                dist += diff / sum;
                        }     
                }
        } else {
                for (int i = 0; i < P_len; i++){
                        diff = fabs(P[i] - Q[i]);
                        sum  = P[i] + Q[i];
                        // replace 0/0 or x/0 or 0/x by 0 according to Sung-Hyuk Cha (2007)
                        if ((diff == 0.0) || (sum == 0.0)){
                                dist += 0.0;
                        } else {
                                dist += diff / sum;
                        }     
                }
        }
        return dist;
}


//' @title Lorentzian distance (lowlevel function)
//' @description The lowlevel function for computing the lorentzian distance.
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
double lorentzian(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA, const Rcpp::String unit){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0.0;
        double diff  = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if (testNA){
                for (int i = 0; i < P_len; i++){
                        if ((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        diff = fabs(P[i] - Q[i]);
                        if (unit == "log"){
                                dist += log(1.0 + diff);
                        } 
                        else if (unit == "log2"){
                                dist += custom_log2(1.0 + diff);
                        }
                        else if (unit == "log10"){
                                dist += custom_log10(1.0 + diff);
                        } else {
                                Rcpp::stop("Please choose from units: log, log2, or log10.");
                        }
                 }
        } else {
                
                for (int i = 0; i < P_len; i++){
                        diff = fabs(P[i] - Q[i]);
                        if (unit == "log"){
                                dist += log(1.0 + diff);
                        } 
                        else if (unit == "log2"){
                                dist += custom_log2(1.0 + diff);
                        }
                        else if (unit == "log10"){
                                dist += custom_log10(1.0 + diff);
                        } else {
                                Rcpp::stop("Please choose from units: log, log2, or log10.");
                        }
                 }
        }
        return dist;
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
double intersection_dist(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        double min_point  = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if (testNA){
                for (int i = 0; i < P_len; i++){
                        if ((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        if (P[i] <= Q[i]){
                                min_point = P[i];
                        } else {
                                min_point = Q[i];
                        }
                        dist += min_point;
                }   
        } else {
                for (int i = 0; i < P_len; i++){
                        if (P[i] <= Q[i]){
                                min_point = P[i];
                        } else {
                                min_point = Q[i];
                        }
                        dist += min_point;
                } 
        }
        return dist;
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
double wave_hedges(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double diff       = 0.0;
        double dist       = 0.0;
        double max_point  = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if (testNA){
                for (int i = 0; i < P_len; i++){
                        if ((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        diff      = fabs(P[i] - Q[i]);
                        if (P[i] >= Q[i]){
                                max_point = P[i];
                        } else {
                                max_point = Q[i];
                        }
                        if ((diff == 0.0) || (max_point == 0.0)){
                                dist += 0.0;
                        } else {
                                dist += diff / max_point;
                        }     
                }
        } else {
                for (int i = 0; i < P_len; i++){
                        diff      = fabs(P[i] - Q[i]);
                        if (P[i] >= Q[i]){
                                max_point = P[i];
                        } else {
                                max_point = Q[i];
                        }
                        if ((diff == 0.0) || (max_point == 0.0)){
                                dist += 0.0;
                        } else {
                                dist += diff / max_point;
                        }     
                }
        }
        return dist;       
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
double czekanowski(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double diff       = 0.0;
        double sum        = 0.0;
        double dist1      = 0.0;
        double dist2      = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if (testNA){
                for (int i = 0; i < P_len; i++){
                        if ((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        diff  = fabs(P[i] - Q[i]);
                        sum   = P[i] + Q[i];
                        dist1 += diff;
                        dist2 += sum;
                }
        } else {
                for (int i = 0; i < P_len; i++){
                        diff  = fabs(P[i] - Q[i]);
                        sum   = P[i] + Q[i];
                        dist1 += diff;
                        dist2 += sum;
                }
        }
        
        if (dist2 == 0.0) {
                return NAN;
        } else {
                return dist1/dist2;         
        }       
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
double motyka(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double sum        = 0.0;
        double dist1      = 0.0;
        double dist2      = 0.0;
        double min_point  = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                if(testNA){
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                }
                
                sum = P[i] + Q[i];
                
                if (P[i] <= Q[i]){
                        
                        min_point = P[i];
                        
                } else {
                        
                        min_point = Q[i];
                }
               
                dist1 += min_point;
                dist2 += sum;
                
        }
        
        if (dist2 == 0.0) {
                return NAN;
        } else {
                return (1.0 - (dist1 / dist2));        
        } 
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
double tanimoto(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        // Soergel = Tanimoto
        return soergel(P, Q, testNA);
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
double ruzicka(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        // Ruzicka = 1 - Tanimoto = 1 - Soergel
        return (1.0 - soergel(P, Q, testNA));
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
double inner_product(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                if(testNA){
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                }
                
                dist += P[i] * Q[i];
        }
        
        return dist;  
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
double harmonic_mean_dist(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len     = P.size();
        int    Q_len     = Q.size();
        double prod      = 0.0;
        double sum       = 0.0;
        double dist      = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                if(testNA){
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                }
                
                prod = P[i] * Q[i];
                sum  = P[i] + Q[i];
                
                if((prod == 0.0) || (sum == 0.0)){
                        dist += 0.0;
                } else {
                        
                        dist += prod / sum;
                } 
        }
        
        return (2.0 * dist);      
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
double cosine_dist(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len     = P.size();
        int    Q_len     = Q.size();
        double prod      = 0.0;
        double p_square  = 0.0;
        double q_square  = 0.0;
        double dist      = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                if(testNA){
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                }
                
                prod      = P[i] * Q[i];
                p_square  += pow(P[i], 2.0);
                q_square  += pow(Q[i], 2.0);
                dist      += prod;
                
        }
        
        double dist2 = sqrt(p_square) * sqrt(q_square);
        
        if (dist2 == 0.0) {
                return NAN;
        } else {
                return dist / dist2;
        }
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
double kumar_hassebrook(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len     = P.size();
        int    Q_len     = Q.size();
        double prod      = 0.0;
        double p_square  = 0.0;
        double q_square  = 0.0;
        double dist      = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                if(testNA){
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                }
                
                prod      = P[i] * Q[i];
                p_square  += pow(P[i], 2.0);
                q_square  += pow(Q[i], 2.0);
                dist      += prod;
                
        }
        
        double dist2 = p_square + q_square - dist;
        
        if (dist2 == 0.0) {
                return 0;
        } else {
                return dist / dist2;
        }
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
double jaccard(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        return (1.0 - kumar_hassebrook(P,Q, testNA));
        
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
double dice_dist(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len       = P.size();
        int    Q_len       = Q.size();
        double diff_square = 0.0;
        double p_square    = 0.0;
        double q_square    = 0.0;
        double dist        = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                if(testNA){
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                }
                
                diff_square =  pow((P[i] - Q[i]), 2.0);
                p_square    += pow(P[i], 2.0);
                q_square    += pow(Q[i], 2.0);
                dist        += diff_square;
                
        }
        
        double dist2 = p_square + q_square;
        
        if (dist2 == 0.0) {
                return NAN;
        } else {
                return (dist / dist2);         
        }
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
double fidelity(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist      = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                if(testNA){
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                }
                
                dist += sqrt(P[i] * Q[i]);
                
        }
        
        return dist;   
}


//' @title Bhattacharyya distance (lowlevel function)
//' @description The lowlevel function for computing the bhattacharyya distance.
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
//' bhattacharyya(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE, unit = "log2")
//' @export
// [[Rcpp::export]]
double bhattacharyya(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA, const Rcpp::String unit){
        
        if (unit == "log"){
                return -log(fidelity(P,Q, testNA));
        }
        
        else if (unit == "log2"){
                return -custom_log2(fidelity(P,Q, testNA));
        }
        
        else if (unit == "log10"){
                return -custom_log10(fidelity(P,Q, testNA));
        } else {
                Rcpp::stop("Please choose from units: log, log2, or log10.");
                return -1.0;
        }
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
double hellinger(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        return 2.0 * sqrt( 1.0 - fidelity(P,Q, testNA));
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
double matusita(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        return sqrt( 2.0 - ( 2.0 * fidelity(P,Q, testNA)));
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
double squared_chord(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                if(testNA){
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                }
                
                dist += pow(sqrt(P[i]) - sqrt(Q[i]), 2.0);
                
        }
        
        return dist;
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
double squared_euclidean(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if(testNA){
                
                for(int i = 0; i < P_len; i++){
                
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                  
                        dist += pow(P[i] - Q[i], 2.0);
                
                }
        } else {
                
                for(int i = 0; i < P_len; i++){
                
                        dist += pow(P[i] - Q[i], 2.0);
                
                }   
        }
        
        return dist;     
}


//' @title Pearson chi-squared distance (lowlevel function)
//' @description The lowlevel function for computing the pearson_chi_sq distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' pearson_chi_sq(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double pearson_chi_sq(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if(testNA){
                
                for(int i = 0; i < P_len; i++){
                
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                
                         if(Q[i] == 0.0){
                        
                                dist += pow(P[i] - Q[i], 2.0) / 0.00001 ;
                         } else {
                        
                                dist += pow(P[i] - Q[i], 2.0) / Q[i];
                         }  
                }
        } else {
                
                for(int i = 0; i < P_len; i++){
                
                         if(Q[i] == 0.0){
                        
                                dist += pow(P[i] - Q[i], 2.0) / 0.00001 ;
                         } else {
                        
                                dist += pow(P[i] - Q[i], 2.0) / Q[i];
                         }  
                }
                
        }
        
        return dist;
}


//' @title Neyman chi-squared distance (lowlevel function)
//' @description The lowlevel function for computing the neyman_chi_sq distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' neyman_chi_sq(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double neyman_chi_sq(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if(testNA){
                
                for(int i = 0; i < P_len; i++){
                
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        
                        if(P[i] == 0.0){
                                
                                dist += pow(P[i] - Q[i], 2.0) / 0.00001;
                        } else {
                        
                                dist += pow(P[i] - Q[i], 2.0) / P[i];
                        }
                }
        } else {
                
                for(int i = 0; i < P_len; i++){
                
                        if(P[i] == 0.0){
                                
                                dist += pow(P[i] - Q[i], 2.0) / 0.00001;
                        } else {
                        
                                dist += pow(P[i] - Q[i], 2.0) / P[i];
                        }
                }
        }
        
        return dist;
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
double squared_chi_sq(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        double PQdiff     = 0.0;
        double PQsum      = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
       
        for(int i = 0; i < P_len; i++){
                if(testNA){
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                }
                
                PQdiff = pow(P[i] - Q[i], 2.0);
                PQsum  = P[i] + Q[i];
                
                if((PQdiff == 0.0) || (PQsum == 0.0)){
                        
                        dist += 0.0;
                } else {
                        
                        dist += PQdiff / PQsum;
                }                
        }
        
        return dist;   
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
double prob_symm_chi_sq(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        return (2.0 * squared_chi_sq(P,Q, testNA));
        
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
double divergence_sq(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        double PQdiff     = 0.0;
        double PQsum      = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                if(testNA){
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                }
                
                PQdiff = pow(P[i] - Q[i], 2.0);
                PQsum  = pow(P[i] + Q[i], 2.0);
                
                if((PQdiff == 0.0) || (PQsum == 0.0)){
                        
                        dist += 0.0;
                } else {
                        
                        dist += PQdiff / PQsum;
                }
        }
        
        return 2.0 * dist;
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
double clark_sq(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        double PQdiff     = 0.0;
        double PQsum      = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }

        for(int i = 0; i < P_len; i++){
                if(testNA){
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                }
                
                PQdiff = fabs(P[i] - Q[i]);
                PQsum  = P[i] + Q[i];
                
                if((PQdiff == 0.0) || (PQsum == 0.0)){
                        
                        dist += 0.0;
                } else {
                        
                        dist += pow(PQdiff / PQsum, 2.0);
                }
        }
        
        return sqrt(dist);
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
double additive_symm_chi_sq(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        double PQsum      = 0.0;
        double PQprod     = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                if(testNA){
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                }
                
                PQsum  = P[i] + Q[i];
                PQprod = P[i] * Q[i];
                
                if((PQsum == 0.0) || (PQprod == 0.0)){
                        
                        dist += 0.0;
                } else {
                        
                        dist += pow(P[i] - Q[i], 2.0) * (PQsum / PQprod);
                }
        }
        
        return dist;
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
//' @author Hajk-Georg Drost
//' @examples
//' kullback_leibler_distance(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE, unit = "log2")
//' @export
// [[Rcpp::export]]
double kullback_leibler_distance(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA, const Rcpp::String unit){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        double PQratio    = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }

        if(testNA){
                for(int i = 0; i < P_len; i++){
                
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        
                        if((P[i] == 0.0) && (Q[i] == 0.0)){
                                dist += 0.0;
                        } else {
                                
                                if(Q[i] == 0.0){
                                        PQratio = P[i] / 0.00001;
                                } else {
                                
                                        PQratio = P[i] / Q[i];
                                }
                                 
                                if (unit == "log"){
                                        if (PQratio == 0.0) {
                                                dist += 0.0;
                                        } else {
                                                dist += P[i] * log(PQratio); 
                                        }         
                                }
                                
                                else if (unit == "log2"){
                                        if (PQratio == 0.0) {
                                                dist += 0.0;
                                        } else {
                                                dist += P[i] * custom_log2(PQratio); 
                                        }
                                }
                                
                                else if (unit == "log10"){
                                        if (PQratio == 0.0) {
                                                dist += 0.0;
                                        } else {
                                                dist += P[i] * custom_log10(PQratio); 
                                        }
                                } else {
                                        Rcpp::stop("Please choose from units: log, log2, or log10.");
                                }
                        }  
                }
        } else {
                
                for(int i = 0; i < P_len; i++){
                        
                        if((P[i] == 0.0) && (Q[i] == 0.0)){
                                dist += 0.0;
                        } else {
                                
                                if(Q[i] == 0.0){
                                        PQratio = P[i] / 0.00001;
                                } else {
                                
                                        PQratio = P[i] / Q[i];
                                }
                        
                                if (unit == "log"){
                                   if (PQratio == 0.0) {
                                           dist += 0.0;
                                   } else {
                                           dist += P[i] * log(PQratio); 
                                   }      
                                }
                                
                                else if (unit == "log2"){
                                        if (PQratio == 0.0) {
                                                dist += 0.0;
                                        } else {
                                                dist += P[i] * custom_log2(PQratio); 
                                        }
                                }
                                
                                else if (unit == "log10"){
                                        if (PQratio == 0.0) {
                                                dist += 0.0;
                                        } else {
                                                dist += P[i] * custom_log10(PQratio); 
                                        }
                                } else {
                                        Rcpp::stop("Please choose from units: log, log2, or log10.");
                                }
                        }  
                }  
        }
        
        return dist;
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
//' @author Hajk-Georg Drost
//' @examples
//' jeffreys(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE, unit = "log2")
//' @export
// [[Rcpp::export]]
double jeffreys(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA, const Rcpp::String unit){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double PQrate     = 0.0;
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if(testNA){
                
                for(int i = 0; i < P_len; i++){
                
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        
                        if(Q[i] == 0.0){
                                PQrate = P[i] / 0.00001;
                        } else {
                                PQrate = P[i] / Q[i];
                        }
                
                        if(PQrate == 0.0){
                                if (unit == "log"){
                                        dist += (P[i] - Q[i]) * log(0.00001);
                                }
                                
                                else if (unit == "log2"){
                                        dist += (P[i] - Q[i]) * custom_log2(0.00001);
                                } 
                                
                                else if (unit == "log10"){
                                        dist += (P[i] - Q[i]) * custom_log10(0.00001);
                                } else {
                                        Rcpp::stop("Please choose from units: log, log2, or log10.");
                                }
                                
                        } else {
                                if (unit == "log"){
                                        dist += (P[i] - Q[i]) * log(PQrate);
                                }
                                
                                else if (unit == "log2"){
                                        dist += (P[i] - Q[i]) * custom_log2(PQrate);
                                } 
                                
                                else if (unit == "log10"){
                                        dist += (P[i] - Q[i]) * custom_log10(PQrate);
                                } else {
                                        Rcpp::stop("Please choose from units: log, log2, or log10.");
                                }
                        }     
                }
        } else {
                
                for(int i = 0; i < P_len; i++){
                        
                        if(Q[i] == 0.0){
                                PQrate = P[i] / 0.00001;
                        } else {
                                PQrate = P[i] / Q[i];
                        }
                
                        if(PQrate == 0.0){
                                if (unit == "log"){
                                        dist += (P[i] - Q[i]) * log(0.00001);
                                }
                                
                                else if (unit == "log2"){
                                        dist += (P[i] - Q[i]) * custom_log2(0.00001);
                                } 
                                
                                else if (unit == "log10"){
                                        dist += (P[i] - Q[i]) * custom_log10(0.00001);
                                } else {
                                        Rcpp::stop("Please choose from units: log, log2, or log10.");
                                }
                        } else {
                        
                                if (unit == "log"){
                                        dist += (P[i] - Q[i]) * log(PQrate);
                                }
                                
                                else if (unit == "log2"){
                                        dist += (P[i] - Q[i]) * custom_log2(PQrate);
                                } 
                                
                                else if (unit == "log10"){
                                        dist += (P[i] - Q[i]) * custom_log10(PQrate);
                                } else {
                                        Rcpp::stop("Please choose from units: log, log2, or log10.");
                                }
                        }     
                }
                
        }
        
        return dist;
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
double k_divergence(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA, const Rcpp::String unit){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if(testNA){
                
                for(int i = 0; i < P_len; i++){
                        if ((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        
                        if ((P[i] == 0.0) && (Q[i] == 0.0)){
                                dist += 0.0;
                        } else {
                                
                                if (unit == "log"){      
                                    dist += (P[i] * log((2.0 * P[i]) / (P[i] + Q[i])));    
                                }
                        
                                else if (unit == "log2"){
                                        dist += (P[i] * custom_log2((2.0 * P[i]) / (P[i] + Q[i])));
                                }
                                
                                else if (unit == "log10"){
                                        dist += (P[i] * custom_log10((2.0 * P[i]) / (P[i] + Q[i])));
                                } else {
                                        Rcpp::stop("Please choose from units: log, log2, or log10.");
                                }
                        }
                }
        } else {
                
                for (int i = 0; i < P_len; i++){
                        
                        if ((P[i] == 0.0) && (Q[i] == 0.0)){
                                dist += 0.0;
                        } else {
                        
                                if (unit == "log"){      
                                    dist += (P[i] * log((2.0 * P[i]) / (P[i] + Q[i])));    
                                }
                        
                                else if (unit == "log2"){
                                        dist += (P[i] * custom_log2((2.0 * P[i]) / (P[i] + Q[i])));
                                }
                                
                                else if (unit == "log10"){
                                        dist += (P[i] * custom_log10((2.0 * P[i]) / (P[i] + Q[i])));
                                } else {
                                        Rcpp::stop("Please choose from units: log, log2, or log10.");
                                }
                         }
                }
        }
        
        return dist;
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
double topsoe(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA, const Rcpp::String unit){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        double PQsum      = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if(testNA){
                for(int i = 0; i < P_len; i++){
                
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        
                        if((P[i] == 0.0) && (Q[i] == 0.0)){
                                dist += 0.0;
                        } else {
                        
                               PQsum = P[i] + Q[i];
                               
                               if (unit == "log"){
                                       dist += ((P[i] * log((2.0 * P[i]) / PQsum )) + (Q[i] * log((2.0 * Q[i]) / PQsum ))); 
                               }
                               
                               else if (unit == "log2"){
                                       dist += ((P[i] * custom_log2((2.0 * P[i]) / PQsum )) + (Q[i] * custom_log2((2.0 * Q[i]) / PQsum ))); 
                               }
                               
                               else if (unit == "log10"){
                                       dist += ((P[i] * custom_log10((2.0 * P[i]) / PQsum )) + (Q[i] * custom_log10((2.0 * Q[i]) / PQsum ))); 
                               } else {
                                       Rcpp::stop("Please choose from units: log, log2, or log10.");
                               }
                        }
                }
        } else {
                
                for(int i = 0; i < P_len; i++){
                
                        if((P[i] == 0.0) && (Q[i] == 0.0)){
                                dist += 0.0;
                        } else {
                        
                               PQsum = P[i] + Q[i];
                               
                               if (unit == "log"){
                                       dist += ((P[i] * log((2.0 * P[i]) / PQsum )) + (Q[i] * log((2.0 * Q[i]) / PQsum ))); 
                               }
                               
                               else if (unit == "log2"){
                                       dist += ((P[i] * custom_log2((2.0 * P[i]) / PQsum )) + (Q[i] * custom_log2((2.0 * Q[i]) / PQsum ))); 
                               }
                               
                               else if (unit == "log10"){
                                       dist += ((P[i] * custom_log10((2.0 * P[i]) / PQsum )) + (Q[i] * custom_log10((2.0 * Q[i]) / PQsum ))); 
                               } else {
                                       Rcpp::stop("Please choose from units: log, log2, or log10.");
                               }
                        }
                }
        }
        
        return dist;
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
double jensen_shannon(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA, const Rcpp::String unit){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double sum1       = 0.0;
        double sum2       = 0.0;
        double PQsum      = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
       
       if(testNA){
               for(int i = 0; i < P_len; i++){
                       if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                               PQsum =   P[i] + Q[i];
                                
                                if (unit == "log"){
                                  if (P[i] == 0.0 || PQsum == 0.0) {
                                    sum1  += 0.0;
                                  } else {
                                    sum1  +=  P[i] * log((2.0 * P[i]) / PQsum);
                                  }
                                  if (Q[i] == 0.0 || PQsum == 0.0) {
                                    sum2  += 0.0;
                                  } else {
                                    sum2  +=  Q[i] * log((2.0 * Q[i]) / PQsum);
                                  }
                                }
                                
                                else if (unit == "log2"){
                                  if (P[i] == 0.0 || PQsum == 0.0) {
                                    sum1  += 0.0;
                                  } else {
                                    sum1  +=  P[i] * custom_log2((2.0 * P[i]) / PQsum);
                                  }
                                  if (Q[i] == 0.0 || PQsum == 0.0) {
                                    sum2  += 0.0;
                                  } else {
                                    sum2  +=  Q[i] * custom_log2((2.0 * Q[i]) / PQsum);
                                  }
                                }
                                
                                else if (unit == "log10"){
                                  if (P[i] == 0.0 || PQsum == 0.0) {
                                    sum1  += 0.0;
                                  } else {
                                    sum1  +=  P[i] * custom_log10((2.0 * P[i]) / PQsum);
                                  }
                                  if (Q[i] == 0.0 || PQsum == 0.0) {
                                    sum2  += 0.0;
                                  } else {
                                    sum2  +=  Q[i] * custom_log10((2.0 * Q[i]) / PQsum);
                                  }
                                } else {
                                        Rcpp::stop("Please choose from units: log, log2, or log10.");
                                }
                        }
        } else {
                
                for(int i = 0; i < P_len; i++){
                       
                                PQsum =   P[i] + Q[i];
                                
                                if (unit == "log"){
                                  if (P[i] == 0.0 || PQsum == 0.0) {
                                    sum1  += 0.0;
                                  } else {
                                    sum1  +=  P[i] * log((2.0 * P[i]) / PQsum);
                                  }
                                  if (Q[i] == 0.0 || PQsum == 0.0) {
                                    sum2  += 0.0;
                                  } else {
                                    sum2  +=  Q[i] * log((2.0 * Q[i]) / PQsum);
                                  }
                                }
                                
                                else if (unit == "log2"){
                                  if (P[i] == 0.0 || PQsum == 0.0) {
                                    sum1  += 0.0;
                                  } else {
                                    sum1  +=  P[i] * custom_log2((2.0 * P[i]) / PQsum);
                                  }
                                  if (Q[i] == 0.0 || PQsum == 0.0) {
                                    sum2  += 0.0;
                                  } else {
                                    sum2  +=  Q[i] * custom_log2((2.0 * Q[i]) / PQsum);
                                  }
                                }
                                
                                else if (unit == "log10"){
                                  if (P[i] == 0.0 || PQsum == 0.0) {
                                    sum1  += 0.0;
                                  } else {
                                    sum1  +=  P[i] * custom_log10((2.0 * P[i]) / PQsum);
                                  }
                                  if (Q[i] == 0.0 || PQsum == 0.0) {
                                    sum2  += 0.0;
                                  } else {
                                    sum2  +=  Q[i] * custom_log10((2.0 * Q[i]) / PQsum);
                                  }
                                } else {
                                        Rcpp::stop("Please choose from units: log, log2, or log10.");
                                }
                        }
                }
              
        return 0.5 * (sum1 + sum2);
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
double jensen_difference(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA, const Rcpp::String unit){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        double PQsum      = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if(testNA){
                
                for(int i = 0; i < P_len; i++){
                
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }

                                        PQsum = P[i] + Q[i];
                                        
                                        if (unit == "log"){
                                          if (PQsum == 0.0 || P[i] == 0.0 || Q[i] == 0.0) {
                                            if (PQsum == 0.0 && P[i] == 0.0 && Q[i] == 0.0){
                                              dist += 0.0 ;
                                            } 
                                            if (P[i] == 0.0 && Q[i] > 0.0 && PQsum > 0.0) {
                                              dist += ((0.0 + (Q[i] * log(Q[i]))) / 2.0) - ((PQsum / 2.0) * log(PQsum / 2.0)) ;
                                            }
                                            
                                            if (P[i] > 0.0 && Q[i] == 0.0 && PQsum > 0.0) {
                                              dist += (((P[i] * log(P[i])) + 0.0 ) / 2.0) - ((PQsum / 2.0) * log(PQsum / 2.0)) ;
                                            }
                                            
                                            if (P[i] == 0.0 && Q[i] == 0.0 && PQsum > 0.0) {
                                              dist += 0.0 - ((PQsum / 2.0) * log(PQsum / 2.0)) ;
                                            }
                                            
                                            if (P[i] > 0.0 && Q[i] > 0.0 && PQsum == 0.0) {
                                              dist += (((P[i] * log(P[i])) + (Q[i] * log(Q[i]))) / 2.0) - 0.0 ;
                                            }
                                            
                                            if (P[i] > 0.0 && Q[i] == 0.0 && PQsum == 0.0) {
                                              dist += (((P[i] * log(P[i])) + 0.0) / 2.0) - 0.0 ;
                                            }
                                            
                                            if (P[i] == 0.0 && Q[i] > 0.0 && PQsum == 0.0) {
                                              dist += ((0.0 + (Q[i] * log(Q[i]))) / 2.0) - 0.0 ;
                                            }
                                            
                                          } else {
                                            dist += (((P[i] * log(P[i])) + (Q[i] * log(Q[i]))) / 2.0) - ((PQsum / 2.0) * log(PQsum / 2.0)) ;
                                          }    
                                                
                                        }
                                        
                                        else if (unit == "log2"){
                                          
                                          if (PQsum == 0.0 || P[i] == 0.0 || Q[i] == 0.0) {
                                            if (PQsum == 0.0 && P[i] == 0.0 && Q[i] == 0.0){
                                            dist += 0.0 ;
                                          } 
                                           if (P[i] == 0.0 && Q[i] > 0.0 && PQsum > 0.0) {
                                             dist += ((0.0 + (Q[i] * custom_log2(Q[i]))) / 2.0) - ((PQsum / 2.0) * custom_log2(PQsum / 2.0)) ;
                                           }
                                           
                                           if (P[i] > 0.0 && Q[i] == 0.0 && PQsum > 0.0) {
                                             dist += (((P[i] * custom_log2(P[i])) + 0.0 ) / 2.0) - ((PQsum / 2.0) * custom_log2(PQsum / 2.0)) ;
                                           }
                                           
                                           if (P[i] == 0.0 && Q[i] == 0.0 && PQsum > 0.0) {
                                             dist += 0.0 - ((PQsum / 2.0) * custom_log2(PQsum / 2.0)) ;
                                           }
                                           
                                           if (P[i] > 0.0 && Q[i] > 0.0 && PQsum == 0.0) {
                                             dist += (((P[i] * custom_log2(P[i])) + (Q[i] * custom_log2(Q[i]))) / 2.0) - 0.0 ;
                                           }
                                           
                                           if (P[i] > 0.0 && Q[i] == 0.0 && PQsum == 0.0) {
                                             dist += (((P[i] * custom_log2(P[i])) + 0.0) / 2.0) - 0.0 ;
                                           }
                                           
                                           if (P[i] == 0.0 && Q[i] > 0.0 && PQsum == 0.0) {
                                             dist += ((0.0 + (Q[i] * custom_log2(Q[i]))) / 2.0) - 0.0 ;
                                           }
                                           
                                          } else {
                                            dist += (((P[i] * custom_log2(P[i])) + (Q[i] * custom_log2(Q[i]))) / 2.0) - ((PQsum / 2.0) * custom_log2(PQsum / 2.0)) ;
                                          }                                        
                                          
                                      }
                                        
                                        else if (unit == "log10"){
                                          
                                          if (PQsum == 0.0 || P[i] == 0.0 || Q[i] == 0.0) {
                                            if (PQsum == 0.0 && P[i] == 0.0 && Q[i] == 0.0){
                                              dist += 0.0 ;
                                            } 
                                            if (P[i] == 0.0 && Q[i] > 0.0 && PQsum > 0.0) {
                                              dist += ((0.0 + (Q[i] * custom_log10(Q[i]))) / 2.0) - ((PQsum / 2.0) * custom_log10(PQsum / 2.0)) ;
                                            }
                                            
                                            if (P[i] > 0.0 && Q[i] == 0.0 && PQsum > 0.0) {
                                              dist += (((P[i] * custom_log10(P[i])) + 0.0 ) / 2.0) - ((PQsum / 2.0) * custom_log10(PQsum / 2.0)) ;
                                            }
                                            
                                            if (P[i] == 0.0 && Q[i] == 0.0 && PQsum > 0.0) {
                                              dist += 0.0 - ((PQsum / 2.0) * custom_log10(PQsum / 2.0)) ;
                                            }
                                            
                                            if (P[i] > 0.0 && Q[i] > 0.0 && PQsum == 0.0) {
                                              dist += (((P[i] * custom_log10(P[i])) + (Q[i] * custom_log10(Q[i]))) / 2.0) - 0.0 ;
                                            }
                                            
                                            if (P[i] > 0.0 && Q[i] == 0.0 && PQsum == 0.0) {
                                              dist += (((P[i] * custom_log10(P[i])) + 0.0) / 2.0) - 0.0 ;
                                            }
                                            
                                            if (P[i] == 0.0 && Q[i] > 0.0 && PQsum == 0.0) {
                                              dist += ((0.0 + (Q[i] * custom_log10(Q[i]))) / 2.0) - 0.0 ;
                                            }
                                            
                                          } else {
                                            dist += (((P[i] * custom_log10(P[i])) + (Q[i] * custom_log10(Q[i]))) / 2.0) - ((PQsum / 2.0) * custom_log10(PQsum / 2.0)) ;
                                          }     
                                        } else {
                                                Rcpp::stop("Please choose from units: log, log2, or log10.");
                                        }
                }
        } else {
                
                for(int i = 0; i < P_len; i++){
                        
                        PQsum = P[i] + Q[i];
                  
                  if (unit == "log"){
                    if (PQsum == 0.0 || P[i] == 0.0 || Q[i] == 0.0) {
                      if (PQsum == 0.0 && P[i] == 0.0 && Q[i] == 0.0){
                        dist += 0.0 ;
                      } 
                      if (P[i] == 0.0 && Q[i] > 0.0 && PQsum > 0.0) {
                        dist += ((0.0 + (Q[i] * log(Q[i]))) / 2.0) - ((PQsum / 2.0) * log(PQsum / 2.0)) ;
                      }
                      
                      if (P[i] > 0.0 && Q[i] == 0.0 && PQsum > 0.0) {
                        dist += (((P[i] * log(P[i])) + 0.0 ) / 2.0) - ((PQsum / 2.0) * log(PQsum / 2.0)) ;
                      }
                      
                      if (P[i] == 0.0 && Q[i] == 0.0 && PQsum > 0.0) {
                        dist += 0.0 - ((PQsum / 2.0) * log(PQsum / 2.0)) ;
                      }
                      
                      if (P[i] > 0.0 && Q[i] > 0.0 && PQsum == 0.0) {
                        dist += (((P[i] * log(P[i])) + (Q[i] * log(Q[i]))) / 2.0) - 0.0 ;
                      }
                      
                      if (P[i] > 0.0 && Q[i] == 0.0 && PQsum == 0.0) {
                        dist += (((P[i] * log(P[i])) + 0.0) / 2.0) - 0.0 ;
                      }
                      
                      if (P[i] == 0.0 && Q[i] > 0.0 && PQsum == 0.0) {
                        dist += ((0.0 + (Q[i] * log(Q[i]))) / 2.0) - 0.0 ;
                      }
                      
                    } else {
                      dist += (((P[i] * log(P[i])) + (Q[i] * log(Q[i]))) / 2.0) - ((PQsum / 2.0) * log(PQsum / 2.0)) ;
                    }    
                    
                  }
                  
                  else if (unit == "log2"){
                    
                    if (PQsum == 0.0 || P[i] == 0.0 || Q[i] == 0.0) {
                      if (PQsum == 0.0 && P[i] == 0.0 && Q[i] == 0.0){
                        dist += 0.0 ;
                      } 
                      if (P[i] == 0.0 && Q[i] > 0.0 && PQsum > 0.0) {
                        dist += ((0.0 + (Q[i] * custom_log2(Q[i]))) / 2.0) - ((PQsum / 2.0) * custom_log2(PQsum / 2.0)) ;
                      }
                      
                      if (P[i] > 0.0 && Q[i] == 0.0 && PQsum > 0.0) {
                        dist += (((P[i] * custom_log2(P[i])) + 0.0 ) / 2.0) - ((PQsum / 2.0) * custom_log2(PQsum / 2.0)) ;
                      }
                      
                      if (P[i] == 0.0 && Q[i] == 0.0 && PQsum > 0.0) {
                        dist += 0.0 - ((PQsum / 2.0) * custom_log2(PQsum / 2.0)) ;
                      }
                      
                      if (P[i] > 0.0 && Q[i] > 0.0 && PQsum == 0.0) {
                        dist += (((P[i] * custom_log2(P[i])) + (Q[i] * custom_log2(Q[i]))) / 2.0) - 0.0 ;
                      }
                      
                      if (P[i] > 0.0 && Q[i] == 0.0 && PQsum == 0.0) {
                        dist += (((P[i] * custom_log2(P[i])) + 0.0) / 2.0) - 0.0 ;
                      }
                      
                      if (P[i] == 0.0 && Q[i] > 0.0 && PQsum == 0.0) {
                        dist += ((0.0 + (Q[i] * custom_log2(Q[i]))) / 2.0) - 0.0 ;
                      }
                      
                    } else {
                      dist += (((P[i] * custom_log2(P[i])) + (Q[i] * custom_log2(Q[i]))) / 2.0) - ((PQsum / 2.0) * custom_log2(PQsum / 2.0)) ;
                    }                                        
                    
                  }
                  
                  else if (unit == "log10"){
                    
                    if (PQsum == 0.0 || P[i] == 0.0 || Q[i] == 0.0) {
                      if (PQsum == 0.0 && P[i] == 0.0 && Q[i] == 0.0){
                        dist += 0.0 ;
                      } 
                      if (P[i] == 0.0 && Q[i] > 0.0 && PQsum > 0.0) {
                        dist += ((0.0 + (Q[i] * custom_log10(Q[i]))) / 2.0) - ((PQsum / 2.0) * custom_log10(PQsum / 2.0)) ;
                      }
                      
                      if (P[i] > 0.0 && Q[i] == 0.0 && PQsum > 0.0) {
                        dist += (((P[i] * custom_log10(P[i])) + 0.0 ) / 2.0) - ((PQsum / 2.0) * custom_log10(PQsum / 2.0)) ;
                      }
                      
                      if (P[i] == 0.0 && Q[i] == 0.0 && PQsum > 0.0) {
                        dist += 0.0 - ((PQsum / 2.0) * custom_log10(PQsum / 2.0)) ;
                      }
                      
                      if (P[i] > 0.0 && Q[i] > 0.0 && PQsum == 0.0) {
                        dist += (((P[i] * custom_log10(P[i])) + (Q[i] * custom_log10(Q[i]))) / 2.0) - 0.0 ;
                      }
                      
                      if (P[i] > 0.0 && Q[i] == 0.0 && PQsum == 0.0) {
                        dist += (((P[i] * custom_log10(P[i])) + 0.0) / 2.0) - 0.0 ;
                      }
                      
                      if (P[i] == 0.0 && Q[i] > 0.0 && PQsum == 0.0) {
                        dist += ((0.0 + (Q[i] * custom_log10(Q[i]))) / 2.0) - 0.0 ;
                      }
                      
                    } else {
                      dist += (((P[i] * custom_log10(P[i])) + (Q[i] * custom_log10(Q[i]))) / 2.0) - ((PQsum / 2.0) * custom_log10(PQsum / 2.0)) ;
                    }     
                  } else {
                    Rcpp::stop("Please choose from units: log, log2, or log10.");
                  }
                }
                
        }
        
        return dist;
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
//' @author Hajk-Georg Drost
//' @examples
//' taneja(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE, unit = "log2")
//' @export
// [[Rcpp::export]]
double taneja(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA, const Rcpp::String unit){
        
        int    P_len       = P.size();
        int    Q_len       = Q.size();
        double dist        = 0.0;
        double PQsum       = 0.0;
        double denominator = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        if(testNA){
                
                for(int i = 0; i < P_len; i++){
                
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        
                        if((P[i] == 0.0) && (Q[i] == 0.0)){
                                dist += 0.0;
                        } else {
                                PQsum = P[i] + Q[i];
                                denominator = (2.0 * sqrt(P[i] * Q[i]));
                                
                                if(denominator == 0.0){
                                        
                                        if (unit == "log"){
                                                dist += (PQsum / 2.0) * log(PQsum / 0.00001);
                                        }
                                        
                                        else if (unit == "log2"){
                                                dist += (PQsum / 2.0) * custom_log2(PQsum / 0.00001);
                                        }
                                        
                                        else if (unit == "log10"){
                                                dist += (PQsum / 2.0) * custom_log10(PQsum / 0.00001);
                                        } else {
                                                Rcpp::stop("Please choose from units: log, log2, or log10.");
                                        }
                                        
                                } else {
                                        
                                        if (unit == "log"){
                                                dist += (PQsum / 2.0) * log(PQsum / denominator);
                                        }
                                        
                                        else if (unit == "log2"){
                                                dist += (PQsum / 2.0) * custom_log2(PQsum / denominator);
                                        }
                                        
                                        else if (unit == "log10"){
                                                dist += (PQsum / 2.0) * custom_log10(PQsum / denominator);
                                        } else {
                                                Rcpp::stop("Please choose from units: log, log2, or log10.");
                                        }
                                }
                        }      
                }
        } else {
                
                for(int i = 0; i < P_len; i++){
                        
                        if((P[i] == 0.0) && (Q[i] == 0.0)){
                                dist += 0.0;
                        } else {
                                
                                PQsum = P[i] + Q[i];
                                denominator = (2.0 * sqrt(P[i] * Q[i]));
                                
                                if(denominator == 0.0){
                                        
                                        if (unit == "log"){
                                                dist += (PQsum / 2.0) * log(PQsum / 0.00001);
                                        }
                                        
                                        else if (unit == "log2"){
                                                dist += (PQsum / 2.0) * custom_log2(PQsum / 0.00001);
                                        }
                                        
                                        else if (unit == "log10"){
                                                dist += (PQsum / 2.0) * custom_log10(PQsum / 0.00001);
                                        } else {
                                                Rcpp::stop("Please choose from units: log, log2, or log10.");
                                        }
                                        
                                } else {
                                        
                                        if (unit == "log"){
                                                dist += (PQsum / 2.0) * log(PQsum / denominator);
                                        }
                                        
                                        else if (unit == "log2"){
                                                dist += (PQsum / 2.0) * custom_log2(PQsum / denominator);
                                        }
                                        
                                        else if (unit == "log10"){
                                                dist += (PQsum / 2.0) * custom_log10(PQsum / denominator);
                                        } else {
                                                Rcpp::stop("Please choose from units: log, log2, or log10.");
                                        }
                                }
                        }      
                }
        }
        
        return dist;
}


//' @title Kumar-Johnson distance (lowlevel function)
//' @description The lowlevel function for computing the kumar_johnson distance.
//' @param P a numeric vector storing the first distribution.
//' @param Q a numeric vector storing the second distribution.
//' @param testNA a logical value indicating whether or not distributions shall be checked for \code{NA} values.
//' @author Hajk-Georg Drost
//' @examples
//' kumar_johnson(P = 1:10/sum(1:10), Q = 20:29/sum(20:29), testNA = FALSE)
//' @export
// [[Rcpp::export]]
double kumar_johnson(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double divisor    = 0.0;
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        
        if(testNA){
                for(int i = 0; i < P_len; i++){
                
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                        
                        divisor = (2.0 * pow(P[i] * Q[i], 1.5));
                        
                        if(divisor == 0.0){
                                dist += pow(pow(P[i], 2.0) - pow(Q[i], 2.0), 2.0) / 0.00001;
                        } else {
                                dist += pow(pow(P[i], 2.0) - pow(Q[i], 2.0), 2.0) / divisor;
                        }
                }
        } else {
                
                for(int i = 0; i < P_len; i++){
                
                        divisor = (2.0 * pow(P[i] * Q[i], 1.5));
                        
                        if(divisor == 0.0){
                                dist += pow(pow(P[i], 2.0) - pow(Q[i], 2.0), 2.0) / 0.00001;
                        } else {
                                dist += pow(pow(P[i], 2.0) - pow(Q[i], 2.0), 2.0) / divisor;
                        }
                }
        }
        
        return dist;
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
double avg(const Rcpp::NumericVector& P, const Rcpp::NumericVector& Q, bool testNA){
        
        int    P_len       = P.size();
        int    Q_len       = Q.size();
        double dist        = 0.0;
        double PQdiff      = 0.0;
        double PQmax       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                if(testNA){
                        if((Rcpp::NumericVector::is_na(P[i])) || (Rcpp::NumericVector::is_na(Q[i]))){
                                Rcpp::stop("Your input vector stores NA values...");
                        }
                }
                
                PQdiff = fabs(P[i] - Q[i]);
                
                if (PQdiff > PQmax)
                     PQmax = PQdiff;
                     
                dist += PQdiff;
                
        }
        
        return (dist + PQmax) / 2.0;
}

#endif // philentropy_Distances_H





