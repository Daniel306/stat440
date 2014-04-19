#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericVector BartlettFactorCpp(int d, int df){
  
  // Initialize A, it starts off with all 0s
  NumericVector A(Dimension(d,d));
  NumericVector Norms = rnorm((d*(d-1)/2));
  int NormsCount = 0;
  // Not sure if calculating a numeric vector of dfs makes it faster
  for(int row = 0; row < d; row++){
    A[row*(d+1)] = sqrt(rchisq(1,d-row)[0]);
    for(int col = 0; col < row; col++){
      A[row*d+col] = Norms[NormsCount++];
    }                   
  }
  
  
  return A; // should probably do it with pointers later
}


// [[Rcpp::export]]
SEXP rNIW_Rcpp_2(int n, int d, NumericVector Mu, int kappa, NumericMatrix gamma_inv, int df) {
  NumericVector V_ans(Dimension(d,d,n));
  NumericVector X_ans(Dimension(d,n));
  for (int k = 0; k < n; k++){
   // First, create A.
  NumericVector A = BartlettFactorCpp(d,df);
  
   for(int i = 0; i < d; i++)
    for(int j = 0; j < d; j++){
     
      
      
      
    }
  }
   
   
   
   List ret;
   ret["V"] = V_ans;
   ret["X"] = X_ans;
   return ret;
}
