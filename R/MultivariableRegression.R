

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
  
  
  mNIW.Mu <- A %*% Lambda + (diag(p)- A) %*% beta.hat;
  mNIW.Psi <- Psi + S + C;
  mNIW.df <- df + n;
  mNIW.Kappa <- Kappa;
  mNIW.n <- samples
  
  result <- rMNIW.Rcpp(mNIW.n, mNIW.Mu, mNIW.Kappa, mNIW.Psi, mNIW.df);
  
}

source("rNIW.R")
source("dNIW.R")
source("testcorrectness.R")
#
#X <- matrix(c(1,2,3,5,8,13),3,2)
#Y <- matrix(c(3,5,8), 3, 1)
#result <- generate.Multivariables(100000,X,Y,matrix(0,1,1),10, matrix(0,2,1),matrix(0,2,2))
#rowMeans(result$X, dims = 2)
#rowMeans(result$V, dims = 2)


mNIW.n <- 10000
mNIW.Mu <- matrix(0, 3,4)
mNIW.Kappa <- diag(3)
mNIW.Psi <- diag(4)
mNIW.df <- 5
result <- rMNIW.Rcpp(mNIW.n, mNIW.Mu, mNIW.Kappa, mNIW.Psi, mNIW.df);

f = IW.marginal(1,1,mNIW.Psi,mNIW.df)
curve(f, from = 0, to =1)
plot.compare(f,result$V[1,1,])
