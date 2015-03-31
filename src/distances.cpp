#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;
using namespace std;

//' @export
// [[Rcpp::export]]
double euclidean(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0;
        double diff  = 0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                diff = P[i] - Q[i];
                dist += diff * diff;
                        
        }
        
        return sqrt(dist);
}


//' @export
// [[Rcpp::export]]
double manhattan(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0;
        double diff  = 0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
          
        for(int i = 0; i < P_len; i++){
                
                diff = fabs(P[i] - Q[i]);
                dist += diff;
                        
        }
        
        return dist;
}



//' @export
// [[Rcpp::export]]
double minkowski(const std::vector<double>& P, const std::vector<double>& Q, double n){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0.0;
        double diff  = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                diff = pow(fabs(P[i] - Q[i]), n);
                dist += diff;
                        
        }
        
        return pow(dist ,  1.0/n);
}


//' @export
// [[Rcpp::export]]
double chebyshev(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0;
        double diff  = 0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                diff = fabs(P[i] - Q[i]);
                
                if (diff > dist)
                     dist = diff;
                        
        }
        
        return dist;
}


//' @export
// [[Rcpp::export]]
double sorensen(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist1 = 0;
        double dist2 = 0;
        double diff  = 0;
        double sum   = 0;
        
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                diff = fabs(P[i] - Q[i]);
                sum  = P[i] + Q[i];
                
                dist1 += diff;
                dist2 += sum;
                
        }
        
        
        return dist1/dist2;
}



//' @export
// [[Rcpp::export]]
double gower(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double diff  = 0;
        double dist  = 0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                diff = fabs(P[i] - Q[i]);
                dist += diff;
                
        }
        
        return ((1.0/P_len) * dist);    
}



//' @export
// [[Rcpp::export]]
double soergel(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double diff       = 0;
        double dist1      = 0;
        double dist2      = 0;
        double max_point  = 0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                diff      = fabs(P[i] - Q[i]);
                
                if (P[i] >= Q[i]){
                        
                        max_point = P[i];
                        
                } else {
                        
                        max_point = Q[i];
                }
               
                dist1 += diff;
                dist2 += max_point;
                
        }
        
        return dist1/dist2;     
}


//' @export
// [[Rcpp::export]]
double kulczynski_d(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double diff       = 0;
        double dist1      = 0;
        double dist2      = 0;
        double min_point  = 0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        
        for(int i = 0; i < P_len; i++){
                
                diff      = fabs(P[i] - Q[i]);
                
                if (P[i] <= Q[i]){
                        
                        min_point = P[i];
                        
                } else {
                        
                        min_point = Q[i];
                }
               
                dist1 += diff;
                dist2 += min_point;
                
        }
        
        return dist1/dist2;      
}


//' @export
// [[Rcpp::export]]
double canberra(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0;
        double diff  = 0;
        double sum   = 0;
        
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
            
        for(int i = 0; i < P_len; i++){
                
                diff = fabs(P[i] - Q[i]);
                sum  = P[i] + Q[i];
                dist += diff / sum;
                
        }
        
        return dist;
}


//' @export
// [[Rcpp::export]]
double lorentzian(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len = P.size();
        int    Q_len = Q.size();
        double dist  = 0;
        double diff  = 0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
            
        for(int i = 0; i < P_len; i++){
                
                diff = fabs(P[i] - Q[i]);
                dist += log(1.0 + diff);
                        
        }
        
        return dist;
}


//' @export
// [[Rcpp::export]]
double intersection_dist(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0;
        double min_point  = 0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                                
                if (P[i] <= Q[i]){
                        
                        min_point = P[i];
                        
                } else {
                        
                        min_point = Q[i];
                }
               
                dist += min_point;
                
        }
        
        return dist;
}


//' @export
// [[Rcpp::export]]
double wave_hedges(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double diff       = 0;
        double dist       = 0;
        double max_point  = 0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                diff      = fabs(P[i] - Q[i]);
                
                if (P[i] >= Q[i]){
                        
                        max_point = P[i];
                        
                } else {
                        
                        max_point = Q[i];
                }
               
                dist += diff / max_point;
                
        }
        
        return dist;       
}


//' @export
// [[Rcpp::export]]
double czekanowski(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double diff       = 0;
        double sum        = 0;
        double dist1      = 0;
        double dist2      = 0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
          
        for(int i = 0; i < P_len; i++){
                
                diff  = fabs(P[i] - Q[i]);
                sum   = P[i] + Q[i];
                dist1 += diff;
                dist2 += sum;
                
        }
        
        return dist1 / dist2;      
}



//' @export
// [[Rcpp::export]]
double motyka(const std::vector<double>& P, const std::vector<double>& Q){
        
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
                
                sum      = P[i] + Q[i];
                
                if (P[i] <= Q[i]){
                        
                        min_point = P[i];
                        
                } else {
                        
                        min_point = Q[i];
                }
               
                dist1 += min_point;
                dist2 += sum;
                
        }
        
        return (1.0 - (dist1 / dist2));  
}


//' @export
// [[Rcpp::export]]
double tanimoto(const std::vector<double>& P, const std::vector<double>& Q){
        
        // Soergel = Tanimoto
        return soergel(P, Q);
}


//' @export
// [[Rcpp::export]]
double ruzicka(const std::vector<double>& P, const std::vector<double>& Q){
        
        // Ruzicka = 1 - Tanimoto = 1 - Soergel
        return (1.0 - soergel(P, Q));
}


//' @export
// [[Rcpp::export]]
double inner_product(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                dist += P[i] * Q[i];
                
        }
        
        return dist;  
}


//' @export
// [[Rcpp::export]]
double harmonic_mean_dist(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len     = P.size();
        int    Q_len     = Q.size();
        double prod      = 0.0;
        double sum       = 0.0;
        double dist      = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                prod = P[i] * Q[i];
                sum  = P[i] + Q[i];
                
                dist += prod / sum;
                
        }
        
        return (2.0 * dist);      
}


//' @export
// [[Rcpp::export]]
double cosine_dist(const std::vector<double>& P, const std::vector<double>& Q){
        
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
                
                prod      = P[i] * Q[i];
                p_square  += pow(P[i], 2);
                q_square  += pow(Q[i], 2);
                dist      += prod;
                
        }
        
        return (dist / (sqrt(p_square) * sqrt(q_square)));  
}


//' @export
// [[Rcpp::export]]
double kumar_hassebrook(const std::vector<double>& P, const std::vector<double>& Q){
        
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
                
                prod      = P[i] * Q[i];
                p_square  += pow(P[i], 2);
                q_square  += pow(Q[i], 2);
                dist      += prod;
                
        }
        
        return (dist / (p_square + q_square - dist));       
}


//' @export
// [[Rcpp::export]]
double jaccard(const std::vector<double>& P, const std::vector<double>& Q){
        
        return (1.0 - kumar_hassebrook(P,Q));
        
}


//' @export
// [[Rcpp::export]]
double dice_dist(const std::vector<double>& P, const std::vector<double>& Q){
        
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
                
                diff_square =  pow((P[i] - Q[i]), 2);
                p_square    += pow(P[i], 2);
                q_square    += pow(Q[i], 2);
                dist        += diff_square;
                
        }
        
        return (dist / (p_square + q_square));   
}


//' @export
// [[Rcpp::export]]
double fidelity(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist      = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                dist += sqrt(P[i] * Q[i]);
                
        }
        
        return dist;   
}


//' @export
// [[Rcpp::export]]
double bhattacharyya(const std::vector<double>& P, const std::vector<double>& Q){
        
        return -log(fidelity(P,Q));
}


//' @export
// [[Rcpp::export]]
double hellinger(const std::vector<double>& P, const std::vector<double>& Q){
        
        return 2.0 * sqrt( 1.0 - fidelity(P,Q));
}

//' @export
// [[Rcpp::export]]
double matusita(const std::vector<double>& P, const std::vector<double>& Q){
        
        return sqrt( 2.0 - ( 2.0 * fidelity(P,Q)));
}

//' @export
// [[Rcpp::export]]
double squared_chord(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                dist += pow(sqrt(P[i]) - sqrt(Q[i]), 2);
                
        }
        
        return dist;
}



//' @export
// [[Rcpp::export]]
double squared_euclidean(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                dist += pow(P[i] - Q[i], 2);
                
        }
        
        return dist;     
}


//' @export
// [[Rcpp::export]]
double pearson_chi_sq(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                dist += pow(P[i] - Q[i], 2) / Q[i];
                
        }
        
        return dist;
}


//' @export
// [[Rcpp::export]]
double neyman_chi_sq(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                dist += pow(P[i] - Q[i], 2) / P[i];
                
        }
        
        return dist;
}


//' @export
// [[Rcpp::export]]
double squared_chi_sq(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                dist += (pow(P[i] - Q[i], 2) / (P[i] + Q[i]));
                
        }
        
        return dist;   
}


//' @export
// [[Rcpp::export]]
double prob_symm_chi_sq(const std::vector<double>& P, const std::vector<double>& Q){
        
        return (2.0 * squared_chi_sq(P,Q));
        
}




//' @export
// [[Rcpp::export]]
double divergence_sq(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                dist += (pow(P[i] - Q[i], 2) / pow(P[i] + Q[i], 2));
                
        }
        
        return 2.0 * dist;
}



//' @export
// [[Rcpp::export]]
double clark_sq(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                dist += pow(fabs(P[i] - Q[i]) / (P[i] + Q[i]), 2);
                
        }
        
        return sqrt(dist);
}


//' @export
// [[Rcpp::export]]
double additive_symm_chi_sq(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                dist += (pow(P[i] - Q[i], 2) * (P[i] + Q[i])) / (P[i] * Q[i]);
                
        }
        
        return dist;
}



//' @export
// [[Rcpp::export]]
double kullback_leibler_distance(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                dist += P[i] * log(P[i] / Q[i]);
                
        }
        
        return dist;
}

//' @export
// [[Rcpp::export]]
double jeffreys(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                dist += (P[i] - Q[i]) * log(P[i] / Q[i]);
                
        }
        
        return dist;
}


//' @export
// [[Rcpp::export]]
double k_divergence(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                dist += (P[i] * log((2.0 * P[i]) / (P[i] + Q[i])));
                
        }
        
        return dist;
}


//' @export
// [[Rcpp::export]]
double topsoe(const std::vector<double>& P, const std::vector<double>& Q){
        
        int    P_len      = P.size();
        int    Q_len      = Q.size();
        double dist       = 0.0;
        double PQsum      = 0.0;
        
        if (P_len != Q_len){
                Rcpp::stop("The vectors you are comparing do not have the same length!");
        }
        
        for(int i = 0; i < P_len; i++){
                
                PQsum = P[i] + Q[i];
                
                dist += ((P[i] * log((2.0 * P[i]) / PQsum )) + (Q[i] * log((2.0 * Q[i]) / PQsum )));
                
        }
        
        return dist;
}
