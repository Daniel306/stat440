
source("rNIW.R")

rMN <- function(n, Mu, U, V) {
  # sample the Matrix Normal distribution
  # XXX this code is partially repeated in rMNIW.Rcpp

  
}

f <- function(X, B, V) {
  # sample from the multivariable normal d
  #
  # the model this generates from is Y = XB + E
  #  with 
  # if q = 1, this is the usual multivariate normal (e.g. mvrnorm())
  #  
  # 
  # args:
  #  X: a c(q, d) matrix
  #  B: the coefficient matrix; a c(d, q) matrix
  #  V: the covariance of elements down the rows(XXX columns?) of Y; this is a c(d,d) matrix.
  #
  #
  # returns:
  #  Y: a c(q, d) matrix
}

ff <- function() {
  # generate multivariable normal test data
  
}

# Can't think of how to name the function
# It generates a q by q by n matrix for Vs and
# p by q by n matrix for the betas.
# it can either generate n or use them, not sure which
# currently only have generating n
generate.Multivariables <- function(n ,X , Y, Psi, df, Lambda, Omega){
  # Multivariable Regression in a Bayes
  # 
  # args:
  #  n: the number of samples to take
  
  q <- dim(Psi)[1];
  p <- dim(Omega)[1];
  
  # TODO: can add some checks here
  
  #X'X is used often, only time X is used by itself is in S
  X.sq <- crossprod(X)
  
  # calculate kappa, also used in calculating C.
  # while C needs inverse, rMNIW already does the inverse, so not doing it here
  Kappa <- X.sq + Omega;

  # compute the matrix
  beta.hat <- backsolve(X.sq, t(X) %*% Y);
  
  temp <- (Y - X %*% beta.hat);
  S <- crossprod(temp)
  
  A <- solve(X.sq + Omega) %*% Omega;
  
  # could swap X'XB with XY, not sure if we should
  C <- t(beta.hat) %*% X.sq %*% beta.hat + t(Lambda) %*% Omega %*% Lambda 
        - t(X.sq %*% beta.hat + Omega %*% Lambda) %*% solve(Kappa) %*% (X.sq %*% beta.hat + Omega %*% Lambda);
  
  # construct the posterior parameters
  mNIW.Mu <- A %*% Lambda + diag(p);
  mNIW.Psi <- Psi + S + C;
  mNIW.df <- df + p; #XXX checkme
  mNIW.Kappa <- Kappa;
  mNIW.n <- samples
  
  return(rMNIW.Rcpp(n, mNIW.Mu, mNIW.Kappa, mNIW.Psi, mNIW.df))
}


test.rMNIW <- function() {
  message("Smoketesting MNIW.Rcpp")
  print(rMNIW.Rcpp(n, matrix(1:6,2,3), diag(2), diag(3), 10))
  message("---------------------------------")
}
#test.rMNIW()


message("Data is:")
message("Y:")
#print(Y)
message("X:")
#print(X)

message("Sampling posterior distribution")
n = 1000
#n = 2 #DEBUG
q = 2 # number of response variates
d = 3 # number of predictor variates
result<- rMNIW.Rcpp(n, matrix(1:(q*d),q,d), diag(q), diag(d), 10) #XXX THIS ISN'T A POSTERIOR; that will come later
#print(result) #DEBUG

message("Estimated B")
# the rMNIW returns (X,V), but in this case, the posterior X *is* B.
B.hat = apply(result$X, -3, mean) #take mean()s across the 3rd dimension;
                          # results in a q x d matrix (after unflattening apply's results)
dim(B.hat) = c(q,d)
print(B.hat)
message() #newline

message("Estimated V")
V.hat = apply(result$V, -3, mean) #take mean()s across the 3rd dimension;
                          # results in a d x d matrix (or it would, but apply flattens its results)
dim(V.hat) = c(d,d)
print(V.hat)
message() #newline

# TODO: compute (with 'quantile') the confidence intervals for each B and V
# additionally, we have a whole sample from which we can do bootstrap-like things, compute functions of the data, etc
# but for this simple test, getting the 