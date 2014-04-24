
source("rNIW.R")

rMN <- function(n, Mu, U, V) {
  # sample the Matrix Normal distribution
  # XXX this code is partially repeated in rMNIW.Rcpp

    
}


require(MASS) #COPYPASTA; TODO: factor
rmvnorm = MASS::mvrnorm #repair naming convention

rmultivariableregression <- function(points, B, V) {
  # simulate from the multivariable normal model.
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
    if(!all(is.finite(m) && m != 0)) { # hack around a bad corner case
      m = rep(1, d)
    } else {
      m = 1/m;
    }
    X = t(rmvnorm(n, rep(0, d), diag(m)))
    dim.stash = dim(X)
    dim(X) = NULL
    X[is.nan(X) | is.infinite(X)] = -pi #deal with more badness
    dim(X) = dim.stash
    message("badness is at")
    print(which(!is.finite(X)))
  } else {
    n = dim(X)[1]
    X = points
  }
  
  E = t(rmvnorm(n, rep(0,q), V))
  if(any(!is.finite(E))) { #introduce bias to quickly repair numerical instability
    dim.stash = dim(E)
    E[!is.finite(E)] = -0.1;
    dim(E) = dim.stash
  }
  Y = X %*% B + E

  list(X = X, Y = Y)
}


lm.multivariable <- function(n ,X , Y, Psi, df, Lambda, Omega) { # TODO: Not sure what to name this.
  # Multivariable Regression via Bayesian Reasoning
  # 
  # This fits the model Y = XB + E where (E_i)^T ~iid MultiNormal(0, V)
  #  but it doesn't fit it in the frequentist sense of constructing
  #  true values for B and V -- instead it produces samples from the
  #  Bayesian Posterior.
  #
  #
  # args:
  #  n: the number of samples to take
  #  X, Y: the observed data; these should be (n,q) and (n,d)
  #  Psi, df, Lambda, Omega: your prior knowledge, encoded
  #    as parameters of the Normal|Inverse-Wishart distribution.
  #  TODO: document how to set them to improper priors
  # 
  # returns:
  #  a list containing
  #   $B
  #   $V
  
  q <- dim(Psi)[1];
  p <- dim(Omega)[1];
  
  # TODO: can add some checks here
  
  #X'X is used often, only time X is used by itself is in S
  X.sq <- crossprod(X)
  
  # calculate kappa, also used in calculating C.
  # while C needs inverse, rMNIW already does the inverse, so not doing it here
  Kappa <- X.sq + Omega;

  # this posterior is neat:
  # it actually involves taking the usual frequentist fit (as produced by lm())
  #  B = (X'X)^{-1}X'y
  # and trading off between that fit and the information from the prior.
  
  
  beta.hat <- backsolve(X.sq, t(X) %*% Y);
  
  residuals <- (Y - X %*% beta.hat);
  S <- crossprod(residuals)
  
  A <- backsolve(X.sq + Omega, Omega);
  
  # could swap X'XB with XY, not sure if we should
  C <- t(beta.hat) %*% X.sq %*% beta.hat + t(Lambda) %*% Omega %*% Lambda 
        - t(X.sq %*% beta.hat + Omega %*% Lambda) %*% solve(Kappa) %*% (X.sq %*% beta.hat + Omega %*% Lambda);
  
  # construct the posterior parameters
  mNIW.Mu <- A %*% Lambda + diag(p);
  mNIW.Psi <- Psi + S + C;
  mNIW.df <- df + p; #XXX checkme
  mNIW.Kappa <- Kappa;
  
  return(rMNIW.Rcpp(n, mNIW.Mu, mNIW.Kappa, mNIW.Psi, mNIW.df))
}


test.rMNIW <- function() {
  message("Smoketesting MNIW.Rcpp")
  print(rMNIW.Rcpp(n, matrix(1:6,2,3), diag(2), diag(3), 10))
  message("---------------------------------")
}
#test.rMNIW()



test.EversonMorris <- function(n=33, m=1000) {
  # Smoketest for lm.multivariable
  #
  # simulate multivariable normal test data
  # then see if the bayesian fitter can pull out the
  # correct (artificial and frequentist) coefficients.
  # Parameters based on Everson & Morris [2002]
  #
  #
  # args:
  #  n: number of observed samples to take
  #     the default is purposely small, to reflect
  #     a realistic data-gathering situation.
  #  m: number of posterior samples to take
  #
  # returns:
  #  nothing; instead, results are printed as work is done.

  d = 9;
  q = 2;
  V = matrix(c(0,0,0,0), q, q)
  B = matrix(rnorm(d*q, 0, 5), d, q)

  message("True B")
  print(B)
  message("True V")
  print(V)
  
  data = rmultivariableregression(n, B, V);

  message("Hiding true values from ourselves")
  rm(B, V, d, q)
  
  message("Data is:")
  message("Y:")
  print(data$Y)
  message("X:")
  print(data$X)

  message("Sampling posterior distribution")
  m = 2 #DEBUG
  result <- lm.multivariable(m, data$X, data$Y, NULL, NULL, NULL, NULL) #FIXME
  
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
  #
  # additionally, we have a whole sample from which we can do bootstrap-like things, compute functions of the data, etc
  # but for this simple test, getting the right coefficients is good enough.
}
#test.EversonMorris()

