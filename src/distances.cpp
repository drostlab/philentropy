#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
double euclidean(NumericVector P, NumericVector Q){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0;
        double diff  = 0;
        
        if(P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        
        for(int i = 0; i < P_len; i++){
                
                diff = P[i] - Q[i];
                dist = dist + (diff * diff);
                        
        }
        
        return sqrt(dist);
}


//' @export
// [[Rcpp::export]]
double manhattan(NumericVector P, NumericVector Q){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0;
        double diff  = 0;
        
        if(P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        
        for(int i = 0; i < P_len; i++){
                
                diff = fabs(P[i] - Q[i]);
                dist = dist + diff;
                        
        }
        
        return dist;
}



//' @export
// [[Rcpp::export]]
double minkowski(NumericVector P, NumericVector Q, double n){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0.0;
        double diff  = 0.0;
        
        if(P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
       
        
        for(int i = 0; i < P_len; i++){
                
                diff = pow(fabs(P[i] - Q[i]), n);
                dist = dist + diff;
                        
        }
        
        return pow(dist ,  1.0/n);
}


//' @export
// [[Rcpp::export]]
double chebyshev(NumericVector P, NumericVector Q){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0;
        double diff  = 0;
        
        if(P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
     
        
        for(int i = 0; i < P_len; i++){
                
                diff = fabs(P[i] - Q[i]);
                
                if(diff > dist)
                     dist = diff;
                        
        }
        
        return dist;
}


//' @export
// [[Rcpp::export]]
double sorensen(NumericVector P, NumericVector Q){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist1 = 0;
        double dist2 = 0;
        double diff  = 0;
        double sum   = 0;
        
        
        if(P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        
        for(int i = 0; i < P_len; i++){
                
                diff = fabs(P[i] - Q[i]);
                sum = P[i] + Q[i];
                
                dist1 = dist1 + diff;
                dist2 = dist2 + sum;
                
        }
        
        
        return dist1/dist2;
}



//' @export
// [[Rcpp::export]]
double gower(NumericVector P, NumericVector Q){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double diff  = 0;
        double dist  = 0;
        
        if(P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        
        for(int i = 0; i < P_len; i++){
                
                diff = fabs(P[i] - Q[i]);
                dist = dist + diff;
                
        }
        
        return ((1.0/P_len) * dist);
        
}



//' @export
// [[Rcpp::export]]
double soergel(NumericVector P, NumericVector Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double diff       = 0;
        double dist1      = 0;
        double dist2      = 0;
        double max_point  = 0;
        
        if(P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        
        for(int i = 0; i < P_len; i++){
                
                diff      = fabs(P[i] - Q[i]);
                
                if(P[i] >= Q[i]){
                        
                        max_point = P[i];
                        
                } else {
                        
                        max_point = Q[i];
                }
               
                dist1 = dist1 + diff;
                dist2 = dist2 + max_point;
                
        }
        
        return dist1/dist2;
        
}


//' @export
// [[Rcpp::export]]
double kulczynski(NumericVector P, NumericVector Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double diff       = 0;
        double dist1      = 0;
        double dist2      = 0;
        double min_point  = 0;
        
        if(P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        
        for(int i = 0; i < P_len; i++){
                
                diff      = fabs(P[i] - Q[i]);
                
                if(P[i] <= Q[i]){
                        
                        min_point = P[i];
                        
                } else {
                        
                        min_point = Q[i];
                }
               
                dist1 = dist1 + diff;
                dist2 = dist2 + min_point;
                
        }
        
        return dist1/dist2;
        
}


//' @export
// [[Rcpp::export]]
double canberra(NumericVector P, NumericVector Q){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist = 0;
        double diff  = 0;
        double sum   = 0;
        
        
        if(P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        
        for(int i = 0; i < P_len; i++){
                
                diff = fabs(P[i] - Q[i]);
                sum  = P[i] + Q[i];
                dist = dist + (diff / sum);
                
        }
        
        
        return dist;
}


//' @export
// [[Rcpp::export]]
double lorentzian(NumericVector P, NumericVector Q){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0;
        double diff  = 0;
        
        if(P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        
        for(int i = 0; i < P_len; i++){
                
                diff = fabs(P[i] - Q[i]);
                
                dist = dist + (log(1.0 + diff));
                        
        }
        
        return dist;
}


//' @export
// [[Rcpp::export]]
double intersection_dist(NumericVector P, NumericVector Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist      = 0;
        double min_point  = 0;
        
        if(P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        
        for(int i = 0; i < P_len; i++){
                                
                if(P[i] <= Q[i]){
                        
                        min_point = P[i];
                        
                } else {
                        
                        min_point = Q[i];
                }
               
                dist = dist + min_point;
                
        }
        
        return dist;
        
}
