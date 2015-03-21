#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace std;

//' @export
// [[Rcpp::export]]
double Ecpp(NumericVector Probabilities){
  int len = Probabilities.size();
  double Entropy = 0;

  for(int i = 0; i < len; i++){
    if(Probabilities[i] > 0){
      Entropy += (Probabilities[i] * (log(Probabilities[i])/log(2)));
    }

    else{
      Entropy += 0;
    }
  }

  return(-Entropy);
}

//' @export
// [[Rcpp::export]]
double JEcpp(NumericVector JointProbabilities){
  int len = JointProbabilities.size();
  double JointEntropy = 0;

     for(int i = 0; i < len; i++){
        if(JointProbabilities[i] > 0){

           JointEntropy += (JointProbabilities[i] * (log(JointProbabilities[i])/log(2)));

        } else{

          JointEntropy += 0.0;

        }
     }

  return(-JointEntropy);
}


//' @export
// [[Rcpp::export]]
double CEcpp(NumericVector JointProbabilities,NumericVector Probabilities){

  double ConditionalEntropy = 0;
  // Using the chain rule: H(X | Y) = H(X,Y) - H(Y)
  // Note: it is important that the Probabilities vector corresponds to Y
  ConditionalEntropy = JEcpp(JointProbabilities) - Ecpp(Probabilities);

  return(ConditionalEntropy);
}


//' @export
// [[Rcpp::export]]
double CrossEntropy(NumericVector P, NumericVector Q) {
     double log_2 = log(2);
     int Psize = P.size();
     int Qsize = Q.size();
     
     // Cross-Entropy = Kullback-Leibler Divergence
     double CE = 0.0;
     if(Psize == Qsize){
         for(int i = 0; i < Psize; i++){
            if((P[i] > 0) && (Q[i] > 0)){
                    
                CE += (P[i] * (log(P[i]/Q[i])/log_2));
                
            } else{
              
                CE += 0.0;
         }
       }
        return CE;
        
     } else{       
             
        return -1;
    }
}

//' @export
// [[Rcpp::export]]
double MIcpp(NumericVector X, NumericVector Y, NumericVector XY){

  double MutualInformation = 0;
  // Using the identity: I(X,Y) = H(X) + H(Y) -  H(X,Y)
  MutualInformation = (Ecpp(X) + Ecpp(Y)) - JEcpp(XY);
  return(MutualInformation);
}


//' @export
// [[Rcpp::export]]
double JensonShannonDivergenceCpp(NumericVector P,NumericVector Q){
  int Psize = P.size();
  int Qsize = Q.size();
  double jsd = -1;

  if(Psize == Qsize){
    NumericVector R(Psize);

    for(int i=0; i < Psize; i++){
      R[i] = (P[i] + Q[i])/2;
    }

    jsd = 0.5 * (CrossEntropy(P,R) + CrossEntropy(Q,R));

  }

  return jsd;
}



//' @export
// [[Rcpp::export]]
NumericMatrix Phylotranscriptomics_Probability(NumericMatrix ExpressionSet, NumericVector X, int n_strata) {

      int nCols = ExpressionSet.ncol();
      int nRows = ExpressionSet.nrow();
      NumericMatrix results(n_strata,nCols);
      
        for(int stage = 0; stage < nCols; stage++) {
           for(int ps = 1; ps <= n_strata; ps++){

            double numerator = 0, divisor = 0;

             for(int gene = 0; gene < nRows; gene++) {

                  if(X[gene] == ps){

                      numerator += ExpressionSet(gene, stage);

                  }

                  divisor  += ExpressionSet(gene, stage);
             }

             results(ps-1,stage) = (double) numerator/divisor;
     }
    }
    return results;
}



//' @export
// [[Rcpp::export]]
NumericMatrix Phylotranscriptomics_JointProbability(NumericMatrix ExpressionSet, NumericVector X, NumericVector Y, int n_strata, int nY){

   int nCols = ExpressionSet.ncol();
   int nRows = ExpressionSet.nrow();
   int iterator = 0;
   NumericMatrix resultMatrix(n_strata*nY,nCols);

   for(int stage = 0; stage < nCols; stage++){
     for(int x_it = 1; x_it <= n_strata; x_it++){
       for(int y_it = 1; y_it <= nY; y_it++){

         double numerator = 0, divisor = 0;

         for(int gene = 0; gene < nRows; gene++){

            if((X[gene] == x_it) && (Y[gene] == y_it)){

               numerator += ExpressionSet(gene,stage);

            }

            divisor  += ExpressionSet(gene,stage);
         }
         iterator = (y_it-1) + ((x_it-1)*nY);
         resultMatrix(iterator,stage) = (double) numerator/divisor;
      }
   }
}
return(resultMatrix);
}

//' @export
// [[Rcpp::export]]
NumericVector Phylotranscriptomics_Entropy(NumericMatrix ExpressionSet, NumericVector X, int n_strata) {
     double log_2 = log(2);
     int nCols = ExpressionSet.ncol();
     NumericVector entropy(nCols);
     NumericMatrix P_x = Phylotranscriptomics_Probability(ExpressionSet, X, n_strata);

     for(int stage = 0; stage < nCols; stage++){
          for(int x_it = 0; x_it < n_strata; x_it++){
               if(P_x(x_it,stage) > 0){
                   entropy[stage] += (P_x(x_it,stage) * ((double) log(P_x(x_it,stage))/log_2));
               } else{

                   entropy[stage] += 0;

               }
          }
      entropy[stage] *= -1;

     }

     return entropy;
}


//' @export
// [[Rcpp::export]]
NumericVector Phylotranscriptomics_JointEntropy(NumericMatrix ExpressionSet, NumericVector X, NumericVector Y, int n_strata, int nY) {
   double log_2 = log(2);
   int nCols = ExpressionSet.ncol();
   NumericMatrix P_xy = Phylotranscriptomics_JointProbability(ExpressionSet,X,Y,n_strata,nY);
   int nProbabilities = P_xy.nrow();
   NumericVector joint_entropy(nCols);

   for(int stage = 0; stage < nCols; stage++){
      for(int it = 0; it < nProbabilities; it++){
          if(P_xy(it,stage) > 0){
             joint_entropy[stage] += (P_xy(it,stage) * ((double) log(P_xy(it,stage))/log_2));
          }
          else{
             joint_entropy[stage] += 0;
          }
      }
      joint_entropy[stage] *= -1;
    }

  return joint_entropy;
}


//' @export
// [[Rcpp::export]]
NumericVector Phylotranscriptomics_ConditionalEntropy(NumericMatrix ExpressionSet, NumericVector X, NumericVector Y, int n_strata, int nY){
  // Conditional Entropy: H(X|Y) = H(Y,X) - H(Y)

  NumericVector H_yx = Phylotranscriptomics_JointEntropy(ExpressionSet,Y,X,nY,n_strata);
  NumericVector H_y = Phylotranscriptomics_Entropy(ExpressionSet,X,n_strata);
  int nStage = H_yx.size();
  NumericVector H_x_given_y(nStage);

  for(int stage = 0; stage < nStage; stage++){
      H_x_given_y[stage] = H_yx[stage] - H_y[stage];
    }

  return H_x_given_y;
}


//' @export
// [[Rcpp::export]]
NumericVector Phylotranscriptomics_MutualInformation(NumericMatrix ExpressionSet, NumericVector X, NumericVector Y, int n_strata, int nY) {
   // Mutual Information : I(X,Y) = H(X) + H(Y) - H(X,Y)
   NumericVector H_x = Phylotranscriptomics_Entropy(ExpressionSet,X,n_strata);
   NumericVector H_y = Phylotranscriptomics_Entropy(ExpressionSet,Y,nY);
   NumericVector H_xy = Phylotranscriptomics_JointEntropy(ExpressionSet,X,Y,n_strata,nY);
   int nStages = H_x.size();
   NumericVector mutual_information(nStages);

   for(int stage= 0; stage < nStages; stage++){
       mutual_information[stage] = (H_x[stage] + H_y[stage]) - H_xy[stage];
   }

  return mutual_information;

}



