

# Can't think of how to name the function
# It generates a q by q by n matrix for Vs and
# p by q by n matrix for the betas.
# it can either generate n or use them, not sure which
# currently only have generating n
generate.Multivariables <- function(samples ,X , Y, Psi, df, Lambda, Omega){
  
  q <- dim(Psi)[1];
  p <- dim(Omega)[1];
  n <- dim(X)[1]; # not 100% sure about this
  # can add some checks here
  
  #X'X is used often, only time X is used by itself is in S
  X.sq <- t(X) %*% X;
  
  # calculate kappa, also used in calculating C.
  # while C needs inverse, rMNIW already does the inverse, so not doing it here
  Kappa <- X.sq + Omega;
  
  beta.hat <- solve(X.sq) %*% (t(X) %*% Y);
  
  temp <- (Y - X %*% beta.hat);
  S <- t(temp) %*% temp;
  
  A <- solve(X.sq + Omega) %*% Omega;
  
  # could swap X'XB with XY, not sure if we should
  C <- t(beta.hat) %*% X.sq %*% beta.hat + t(Lambda) %*% Omega %*% Lambda 
        - t(X.sq %*% beta.hat + Omega %*% Lambda) %*% solve(Kappa) %*% (X.sq %*% beta.hat + Omega %*% Lambda);
  
  
  mNIW.Mu <- A %*% Lambda + diag(p);
  mNIW.Psi <- Psi + S + C;
  mNIW.df <- df + n;
  mNIW.Kappa <- Kappa;
  mNIW.n <- samples
  
  result <- rMNIW.Rcpp(mNIW.n, mNIW.Mu, mNIW.Kappa, mNIW.Psi, mNIW.df);
  
}

source("rNIW.R")
source("dNIW.R")
result <- generate.Multivariables(1000,X,Y,diag(2),10, diag(2),diag(2))
rowMeans(result$X, dims = 2)

result<- rMNIW.Rcpp(1000, matrix(1:6,2,3), diag(2), diag(3), 10)
rowMeans(result$X, dims = 2)
