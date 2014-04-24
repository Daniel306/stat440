
source("rNIW.R")

rMN <- function(n, Mu, U, V) {
  # sample the Matrix Normal distribution
  # XXX this code is partially repeated in rMNIW.Rcpp

    
}


require(MASS) #COPYPASTA; TODO: factor
rmvnorm = MASS::mvrnorm #repair naming convention

rmultivariableregression <- function(points, B, V) {
  # sample from the multivariable normal model
  #
  # The model this generates from is Y = XB + E where (E_i)^T ~iid MultiNormal(0, V)
  # which is like the usual least squares model, but now with
  #  multiple columns of Y (and hence, multiple columns of coefficients B)
  # In detail:
  #  let q be the number of response variates
  #  let d be the number of predictors
  #  let n be the number of (i.i.d.) samples
  #  then
  #       Y a c(n, q)
  #       X a c(n, d)
  #       B a c(d, q)
  #       E a c(n, q)
  #         E_i each *row* of E a (1, q)
  #       V is c(q,q) and gives the covariance *within* each single multivariate sample.
  #       
  # if q = 1, this is the usual univariate regression model with normal (e.g. mvrnorm())
  #  
  # 
  # args:
  #  points: either
  #     a c(q, d) matrix X, giving the predictor points to pretend to do measurements at
  #     or, as a convenience,
  #     a scalar n, then number of samples, which is used to construct X
  #           #XXX what about the degenerate case when q=d=1?
  #  B: the coefficient matrix; a c(d, q) matrix
  #  V: the covariance of elements down the rows(XXX columns?) of Y; this is a c(d,d) matrix.
  #
  #
  # returns:
  #  a list with components
  #  X: a c(n, d)
  #  Y: a c(n, q) matrix
  
  if(length(V) == 1) {
    V = matrix(V)
  }
  stopifnot(is.symmetric(V))
  # coerce B t
  if(is.vector(B)) {
     dim(B) = c(length(B), 1)
  }
  d = dim(B)[1]
  q = dim(B)[2]
  
  stopifnot(q == dim(V)[1] && q == dim(V)[2])
  
  
  if(length(points) == 1) {
    n = points
    # special case: generate X sensibly
    # since we don't know where to sample, we might as well sample
    #  i. near the origin
    # ii. such that X is roughly on a scale that evens out
    m = apply(B, 1, mean) #this computes the mean size of the coefficient for each predictor 
    stopifnot(length(m) == d) # DEBUG
    
    X = t(rmvnorm(n, rep(0, d), diag(1/m)))
  } else {
    n = dim(X)[1]
    X = points
  }
  
  E = t(rmvnorm(n, rep(0,q), V))
  Y = X %*% B + E
  return(Y)
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