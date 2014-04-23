#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// // [[Rcpp::export]]
void BartlettFactorCpp(int d, double df, NumericVector &A){
  
  // Initialize A, it starts off with all 0s
  //NumericVector A(Dimension(d,d));
  NumericVector Norms = rnorm((d*(d-1)/2));
  int NormsCount = 0;
  // Not sure if calculating a numeric vector of dfs makes it faster, could do it later.
  for(int col = 0; col < d; col++){
    A[col*(d+1)] = sqrt(rchisq(1,df-col)[0]);
    for(int row = 0; row < col; row++){ 
      A[col*d+row] = Norms[NormsCount++];
    }                   
  }
  
  
  //return A; // should probably do it with pointers later
}


// Used only for inverse on upper triangular matrices
// // [[Rcpp::export]]
NumericVector backSolveInverse(NumericVector &A, int d){
 //
 NumericVector A_inv(Dimension(d,d)); // sets ot to 0

  // preset diagonals in A_inv to 1 to account for subtraction from I
  for(int i = 0; i< d; i++){
    A_inv[i*(d+1)] = 1;
  }


  for(int col = d-1; col > -1; col--){ // colums of solution matrix
    for(int row = d-1; row > -1; row--){ // row of solution matrix, col of A

     for (int row2 = d-1; row2 > row; row2--){ // row of A, row of solution
      A_inv[col*d + row] -= A[row2*d + row]*A_inv[col*d+row2] ;
      }
      A_inv[col*d + row] = A_inv[col*d + row]/A[row*(d+1)]; 
    }
  }
  
  return A_inv;
}


// [[Rcpp::export]]
SEXP rNIW_Rcpp_2(int n, int d, NumericVector Mu, double kappa, NumericVector gamma_inv, double df) {
  NumericVector V_ans(Dimension(d,d,n));
  NumericVector X_ans(Dimension(d,n));
  
  
  NumericVector A(Dimension(d,d));
  for (int k = 0; k < n; k++){
   // First, create A.
  BartlettFactorCpp(d,df, A);
    // Apply backsolve
  NumericVector A_inv = backSolveInverse(A,d);
  

  
  NumericVector z = rnorm(d);
  
  
  NumericVector U_inv(Dimension(d,d));
  // might want to move this to helper function for upper triangular dot product
  for (int col = 0; col < d; col++){ // col of result 
    for(int row = 0; row < d; row++){ // row of result
      for(int i = 0; i < d; i++){
        U_inv[col*d+row] +=  gamma_inv[i*d + row] * A_inv[col*d + i]; // Will optimize after
      }
               // Rcout << U_inv[col*d+row] << " ";
    }  
 //   Rcout << std::endl;
  }
  
  for (int col = 0; col < d; col++){
    for (int row = 0; row < d; row++){
      for(int i = 0; i < d; i++){
        V_ans[k*d*d + col*d + row] += U_inv[i*d + row] * U_inv[i*d + col];  
      }
    }
  }
  
  for(int i = 0; i < d; i++){
    for(int j = 0; j < d; j++){
      X_ans[k*d+i] += U_inv[j*d + i]*z[j];
    }
      X_ans[k*d + i] = X_ans[k*d + i]/sqrt(kappa) + Mu[i];
  }
  
  }
   List ret;
   ret["X"] = X_ans;
   ret["V"] = V_ans;
   return ret;
}
