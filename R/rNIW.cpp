/* rNIW.cpp
 *
 * C-speed routines to back rNIW.R
 *
 * There are two APIs in use here: Rcpp and Eigen.
 * This makes the code more verbose than what its really doing.
 */
#include <Rcpp.h>
using namespace Rcpp;


#include <RcppEigen.h>
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]


// forward declarations
void BartlettFactorCpp(int d, double df, NumericVector &A);


/* References:
 * 
 * - http://stackoverflow.com/questions/15263996/rcpp-how-to-generate-random-multivariate-normal-vector-in-rcpp
 *
 */

/* TODO:
 *
 * - investigate RcppArmadillo <http://dirk.eddelbuettel.com/code/rcpp.armadillo.html> which gives more readable LAPACK calls, among other things.
 * - investigate RcppEigen <http://cran.r-project.org/web/packages/RcppEigen/index.html>
 *   long versus thread on the Rcpp mailing list: http://thread.gmane.org/gmane.comp.lang.r.rcpp/3522
 */


/* from the RcppEigen vignette */
/* this feels like a dirty hack */
MatrixXd crossprod(const MatrixXd& A) {
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>().rankUpdate(A.adjoint());
}

MatrixXd tcrossprod(const MatrixXd& A) {
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>().rankUpdate(A);
}


// [[Rcpp::export]]
NumericMatrix riwishart(NumericMatrix Psi, double df) {

  unsigned int d = Psi.nrow();

  // invert Psi  
  MatrixXd Psi_inv = as<Map<MatrixXd> >(Psi).inverse();
  // TODO: exploit that we know Psi is symmetric (but SelfAdjointView does not contain inverse()!! maybe knowing the matrix is symmetric doesn't help, in that case??)
  // factor Psi
  MatrixXd Gamma = Psi_inv.llt().matrixU();
    
  // step 1: generate the bartlett factor..
  NumericVector A = NumericVector(Dimension(d,d));
  BartlettFactorCpp(d, df, A);

  // step 2: scale it
  MatrixXd A_ = as<Map<MatrixXd> >(A);  //map into Eigen so we can use linear algebras
  MatrixXd U  = A_*Gamma;

  // construct the final symmetric covariance matrix
  //V = (Gamma^T A^T A Gamma)^{-1} = (A Gamma)^{-1} (A Gamma)^{-1}^T = tcrossprod((A Gamma)^{-1})
  MatrixXd V = tcrossprod(U.inverse());
  return Rcpp::wrap(V); 
}

// [[Rcpp::export]]
NumericMatrix riwishart_n(int n, NumericVector Psi, double df) {
  // faster version which computes several Vs at once but only factors Psi once
}


// [[Rcpp::export]]
NumericVector rmvnorm(NumericVector Mu, NumericMatrix V) {
  unsigned int d = Mu.size();
  
  // take the lower chol() of V:
  //   vv^T = V
  // we need lower form because var(v z ) = v var(z) v^T = v I v^T
  MatrixXd v = as<Map<MatrixXd> >(V).llt().matrixL(); //??
  
  // sample standard normals
  NumericVector z = rnorm(d);
  
  //scale the normals to have covariance matrix V
  VectorXd X = v*as<Map<VectorXd> >(z);
  
  // center the normals
  X += as<Map<VectorXd> >(Mu);
  
  return Rcpp::wrap(X);
}


// [[Rcpp::export]]
List rNIW_extremelynaive_eigen(unsigned int n, NumericVector Mu, double kappa, NumericMatrix Psi, double df) {
  // precondition: all the preconditions have been properly checked by rNIW.typechecked
  unsigned int d = Mu.length();
  
  NumericVector V(Dimension(d,d,n));
  NumericVector X(Dimension(d,n));
  
  for(unsigned int i=0; i<n; i++) {
    NumericMatrix Vi = riwishart(Psi, df);

    for(int j=0; j<d; j++) {
    for(int k=0; k<d; k++) {
      V[d*d*i + j*d + k] = Vi(j, k); //XXX this might be the wrong order (but V is symmetric so we can't tell)
    }
    }

    // hack: Rcpp "helpfully" has a whole slew of pseudo-symbolic classes ("Times_Vector_Vector") that it spits out when you try to do operands
    //  I can't tell if you can actually make it collapse things down
    //  So i'm going to use Eigen as a workaround
    // changing the definition of Vi in the process
    Vi = Rcpp::wrap(as<Map<MatrixXd> >(Vi) / kappa);
    //Rcout << "Vi[1,2] = " << Vi[1,2] << std::endl; //DEBUG
    NumericVector Xi = rmvnorm(Mu, Vi);
    for(int j=0; j<d; j++) {
      X[d*i + j] = Xi(j);
    }
    
  }

  return List::create(Named("X") = X, Named("V") = V);
}


// // [[Rcpp::export]] //<-- not exported since A was converted to a reference
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



// now that I can do that.. what?
// [[Rcpp::export]]
List rNIW_snappy_eigen(unsigned int n, NumericVector Mu, double kappa, NumericVector Psi, double df) {
  // precondition: all the preconditions have been properly checked by rNIW.typechecked
  unsigned int d = Mu.length();

  NumericVector V(Dimension(d,d,n));
  NumericVector X(Dimension(d,n));

  // invert Psi  
  MatrixXd Psi_inv = as<Map<MatrixXd> >(Psi).inverse();
  // TODO: exploit that we know Psi is symmetric (but SelfAdjointView does not contain inverse()!! maybe knowing the matrix is symmetric doesn't help, in that case??)
  // factor Psi
  MatrixXd Gamma = Psi_inv.llt().matrixU();

  NumericVector A = NumericVector(Dimension(d,d));
  for(unsigned int i=0; i<n; i++) {
    // step 1: generate the bartlett factor..
    BartlettFactorCpp(d, df, A);
    MatrixXd A_ = as<Map<MatrixXd> >(A);  //map into Eigen so we can use linear algebras
    //TriangularView<MatrixXd, Upper> trA_ = A_.triangularView<Upper>();
    MatrixXd U = A_*Gamma;
    MatrixXd U_inv = U.inverse();
    MatrixXd Vi = tcrossprod(U_inv);
    // memcpy the result out of Eigen and into Rcpp
    // XXX does this need to be done manually? can't we call some "copy this block of memory" operation?
    for(int j=0; j<d; j++) {
    for(int k=0; k<d; k++) {
      V[d*d*i + j*d + k] = Vi(j, k); //XXX this might be the wrong order (but V is symmetric so we can't tell)
    }
    }
    
    // step 2: generate X
    VectorXd Xi = U.triangularView<Upper>().solve(as<Map<VectorXd> >(rnorm(d))); // backsolve(U, z)
    Xi /= sqrt(kappa);
    Xi += as<Map<VectorXd> >(Mu);
    for(int j=0; j<d; j++) {
      X[d*i + j] = Xi(j);
    }
    
  }  

  return List::create(Named("X") = X, Named("V") = V
       //extra debugging stuff
       //, Named("Psi.inv") = Psi_inv
       //, Named("Gamma") = Gamma
       );
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
