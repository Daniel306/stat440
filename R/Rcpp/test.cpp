#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar


// [[Rcpp::export]]
NumericVector Mmult(NumericVector M1, NumericVector M2){
  // R is weird:
  NumericVector dM1 = M1.attr("dim");
  NumericVector dM2 = M2.attr("dim");
  int rowF = dM1[0];
  int colF = dM2[1];
  int d = dM1[1];
  //assert(d = dM2[0]); not sure how assert is
  
  NumericVector ans(Dimension(rowF, colF));
  for (int col = 0; col < colF; col++){
    for (int row = 0; row < rowF; row++){
      for( int i = 0; i < d; i++){
        ans[col*rowF + row] += M1[i*rowF + row] * M2[col*d + i];
      } 
    }
  }
  return ans;
}


// [[Rcpp::export]]
NumericVector generate_z(int d, int p){
  NumericVector norms = rnorm(d*p);
  NumericVector ans(Dimension(d,p)); // might be off here
  //int count = 0;
  for (int col = 0; col < p; col++){
    for(int row = 0; row < d; row++){
      ans[col*d + row] = norms[col*d + row];
    }
  }
  
  return ans;
}
