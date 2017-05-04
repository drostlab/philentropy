//  Part of the philentropy package
//
//  Copyright (C) 2015 - 2017 Hajk-Georg Drost
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

#ifndef philentropy_Correlation_H
#define philentropy_Correlation_H philentropy_Correlation_H

#include <Rcpp.h> 
#include <math.h>
#include <iostream>
#include "utils.h"
#include "distances.h"

// @export
// [[Rcpp::export]]
double pearson_corr_centred(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, bool testNA){
        
        if(x.size() != y.size()){
                Rcpp::stop("Length of input vectors x and y differ!");
        }
        
        int VecLen = x.size();
        // Pearson Correlation Coefficient (Centred)
        double r = 0.0;
        double mean_x = Rcpp::mean(x);
        double mean_y = Rcpp::mean(y);
        Rcpp::NumericVector x_diff (VecLen);
        Rcpp::NumericVector y_diff (VecLen);
        Rcpp::NumericVector x_diff_sq (VecLen);
        Rcpp::NumericVector y_diff_sq (VecLen);
        
        if(testNA){
                for(int i = 0; i < x.size(); i++){
                        
                        if(Rcpp::NumericVector::is_na(x[i]) | Rcpp::NumericVector::is_na(y[i])){
                                Rcpp::stop("Your input vectors store NA values...");
                        }
                        
                        x_diff[i] = x[i] - mean_x;
                        x_diff_sq[i] = pow(x_diff[i], 2.0);
                        y_diff[i] = y[i] - mean_y;
                        y_diff_sq[i] = pow(y_diff[i], 2.0);      
                        }
        } else {
                for(int i = 0; i < x.size(); i++){
                        x_diff[i] = x[i] - mean_x;
                        x_diff_sq[i] = pow(x_diff[i], 2.0);
                        y_diff[i] = y[i] - mean_y;
                        y_diff_sq[i] = pow(y_diff[i], 2.0);
                
                }
        }
                
        r = ( Rcpp::sum(x_diff * y_diff) ) / (sqrt(Rcpp::sum(x_diff_sq)) * sqrt(Rcpp::sum(y_diff_sq)));
        
        return r;
}


// @export
// [[Rcpp::export]]
double pearson_corr_uncentred(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, bool testNA){
        
        if(x.size() != y.size()){
                Rcpp::stop("Length of input vectors x and y differ!");
        }
        
        int VecLen = x.size();
        // Pearson Correlation Coefficient (Uncentred)
        double r_un = 0.0;
        double mean_x = Rcpp::mean(x);
        double mean_y = Rcpp::mean(y);
        Rcpp::NumericVector x_diff (VecLen);
        Rcpp::NumericVector y_diff (VecLen);
        Rcpp::NumericVector x_diff_sq (VecLen);
        Rcpp::NumericVector y_diff_sq (VecLen);
        
        if(testNA){
                for(int i = 0; i < x.size(); i++){
                       if(Rcpp::NumericVector::is_na(x[i]) | Rcpp::NumericVector::is_na(y[i])){
                                Rcpp::stop("Your input vectors store NA values...");
                        } 
                        x_diff[i] = x[i] - mean_x;
                        x_diff_sq[i] = pow(x_diff[i], 2.0);
                        y_diff[i] = y[i] - mean_y;
                        y_diff_sq[i] = pow(y_diff[i], 2.0);
                }
        } else {
                for(int i = 0; i < x.size(); i++){
                        x_diff[i] = x[i] - mean_x;
                        x_diff_sq[i] = pow(x_diff[i], 2.0);
                        y_diff[i] = y[i] - mean_y;
                        y_diff_sq[i] = pow(y_diff[i], 2.0);
                }
        }
        
        r_un = ( Rcpp::sum(Rcpp::NumericVector(x) * Rcpp::NumericVector(y)) ) / (sqrt(Rcpp::sum(x_diff_sq)) * sqrt(Rcpp::sum(y_diff_sq)));
        
        return r_un;
}


// @export
// [[Rcpp::export]]
double squared_pearson_corr(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, bool testNA){
        
        if(x.size() != y.size()){
                Rcpp::stop("Length of input vectors x and y differ!");
        }
        
        int VecLen = x.size();
        // Pearson Correlation Coefficient (Uncentred)
        double r_sq = 0.0;
        double mean_x = Rcpp::mean(x);
        double mean_y = Rcpp::mean(y);
        Rcpp::NumericVector x_diff (VecLen);
        Rcpp::NumericVector y_diff (VecLen);
        Rcpp::NumericVector x_diff_sq (VecLen);
        Rcpp::NumericVector y_diff_sq (VecLen);
        
        if(testNA){
                for(int i = 0; i < x.size(); i++){
                        if(Rcpp::NumericVector::is_na(x[i]) | Rcpp::NumericVector::is_na(y[i])){
                                Rcpp::stop("Your input vectors store NA values...");
                        } 
                        x_diff[i] = x[i] - mean_x;
                        x_diff_sq[i] = pow(x_diff[i], 2.0);
                        y_diff[i] = y[i] - mean_y;
                        y_diff_sq[i] = pow(y_diff[i], 2.0);
                }
        } else {
                for(int i = 0; i < x.size(); i++){
                        x_diff[i] = x[i] - mean_x;
                        x_diff_sq[i] = pow(x_diff[i], 2.0);
                        y_diff[i] = y[i] - mean_y;
                        y_diff_sq[i] = pow(y_diff[i], 2.0); 
                }
        }
        
        r_sq = ( Rcpp::sum(x_diff * y_diff) ) / (sqrt(Rcpp::sum(x * x)) * sqrt(Rcpp::sum(y * y)));
        
        return r_sq * r_sq;
}

//double spearman(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y){
        
//        return pearson_corr_centred(rank(x),rank(y));
//}

#endif // philentropy_Correlation_H

