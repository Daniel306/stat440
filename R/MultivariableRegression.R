
source("rNIW.R")

rMN <- function(n, Mu, U, V) {
  # sample the Matrix Normal distribution
  # XXX this code is partially repeated in rMNIW.Rcpp

    
}


require(MASS) #COPYPASTA; TODO: factor
rmvnorm = MASS::mvrnorm #repair naming convention

rmultivariableregression <- function(points, B, V) { #<-- TODO: rename? r.multivariable.model?
  # simulate from the multivariable normal model.
  #
  # The model this generates from is Y = XB + E where (E_i)^T ~iid MultiNormal(0, V)
  # which is like the usual least squares model, but now with
  #  multiple columns of Y (and hence, multiple columns of coefficients B)
  # In detail:
  #  let q be the number of response variates
  #  let p be the number of predictors
  #  let n be the number of (i.i.d.) samples
  #  then
  #       Y a c(n, q)
  #       X a c(n, p)
  #       B a c(p, q)
  #       E a c(n, q)
  #         E_i each *row* of E a (1, q)
  #       V is c(q,q) and gives the covariance *within* each single multivariate sample.
  #       
  # if q = 1, this is the usual univariate regression model with normal (e.g. mvrnorm())
  #  
  # 
  # args:
  #  points: either
  #     a c(n, p) matrix X, giving the predictor points to pretend to do measurements at
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
    X = rmvnorm(n, rep(0, d), diag(m)) #MASS's mvrnorm returns a n x d matrix, which is just what we need
    dim.stash = dim(X)
    dim(X) = NULL
    X[is.nan(X) | is.infinite(X)] = -pi #deal with more badness
    dim(X) = dim.stash
  } else {
    n = dim(X)[1]
    X = points
  }
  
  stopifnot(dim(X)[1] == n)
  stopifnot(dim(X)[2] == d)
  
  E = rmvnorm(n, rep(0,q), V) #MASS's mvrnorm returns a n x d matrix, which is just what we need
  if(any(!is.finite(E))) { #introduce bias to quickly repair numerical instability
    dim.stash = dim(E)
    E[!is.finite(E)] = -0.1;
    dim(E) = dim.stash
  }
  Y = X %*% B + E

  list(X = X, Y = Y)
}

# TODO: r.random.effects(n, V, W) {
#  # generate samples from the model Y ~ N(mu, V) where itself mu ~ N(


lm.multivariable <- function(m, X, Y, Lambda=NULL, Omega=NULL, Psi=NULL, df=71) {
  # Multivariable Regression via Bayesian Reasoning
  # 
  # This fits the model Y = XB + E where (E_i)^T ~iid MultiNormal(0, V)
  #  but it doesn't fit it in the frequentist sense of constructing
  #  true values for B and V -- instead it produces samples from the
  #  Bayesian Posterior.
  # The posterior, like a good little bayesian, is conjugate in this case:
  #  both the prior and the posterior are Matrix-Normal|Inverse-Wisharts.
  #
  # args:
  #  m: the number of samples to take
  #  X, Y: the observed data; these should be (n,q) and (n,d)
  #  Psi, df, Lambda, Omega: your prior knowledge, encoded
  #    as parameters of the Normal|Inverse-Wishart distribution.
  #   You should provide a df (though this won't make a huge difference)
  #    but, for convenience, by default the other parameters are set
  #    to give "flat" (aka uninformative aka improper) priors. 
  #   The flat prior on Omega causes Lambda to be ignored
  #    (hence its default value: NULL), and, as a special exception,
  #    causes the degrees of freedom of the posterior to change
  #     from (df+n) to (df+n-d).
  #  The flat prior on Psi is orthogonal to the flat prior on Omega (XXX is this true? surely it has some effect, even if only numerical/speed/something)
  #  (the argument order is chosen to reflect the distribution,
  #  which is the result of sampling
  #      MatrixNormal(Lambda, Omega, InverseWishart(Psi, df)), 
  #   though it does makes calling this function awkward)
  # 
  # returns:
  #  a list containing the posterior samples:
  #   $B
  #   $V
  
  n <- dim(X)[1];
  stopifnot(n == dim(Y)[1])
  q <- dim(Y)[2];
  d <- dim(X)[2];
  
  # TODO: can add some checks here
  
  #X'X is used often, only time X is used by itself is in S
  X.sq <- crossprod(X)
  
  # this posterior is neat:
  # it actually involves taking the usual frequentist fit (as produced by lm())
  #  B = (X'X)^{-1}X'y
  # and trading off between that fit and the information from the prior.
  # So, compute the OLS fit beta.hat and its squared-error (SSE) S 
  beta.hat <- solve(X.sq, t(X) %*% Y);
  residuals <- (Y - X %*% beta.hat);
  S <- crossprod(residuals)
  
  # now incorporate the prior parameters
  # the flat priors are special-cased to avoid
  # numerical instability and R's thorny gorgon's type system
  
  if(is.null(Psi)) { Psi = 0 } #let scalar splaying sort it out
  
  if(is.null(Omega)) {
    mNIW.Mu <- beta.hat
    mNIW.Kappa <- X.sq;
    mNIW.Psi <- Psi + S;
    mNIW.df <- df + n - d;
  } else {
    
    if(is.null(Lambda)) {
      stop("You must specify Lambda if you use a non-flat prior on Omega");
    }
    
    # calculate kappa, also used in calculating C.
    # while C needs inverse, rMNIW already does the inverse, so not doing it here
    
    mNIW.Kappa <- X.sq + Omega;
    mNIW.df <- df + n - d;
    
    A <- solve(mNIW.Kappa, Omega);    
    I = diag(d) #identity matrix
    mNIW.Mu <- A %*% Lambda  +  (I-A) %*% beta.hat
    
    # could swap X'XB with XY, not sure if we should
    # this formula is long and grueling
    L = X.sq %*% beta.hat  + Omega %*% Lambda  #this term is used twice
    
    # in all lines below, mahalanobis gives this error 
    # Error in x %*% cov : non-conformable arguments
    C <- t(beta.hat) %*% X.sq %*% beta.hat
           # mahalanobis(beta.hat, 0, X.sq, inverted=TRUE)
             + t(Lambda) %*% Omega %*% Lambda
             #+ mahalanobis(Lambda, 0, Omega, inverted=TRUE) 
             - t(L) %*% solve(mNIW.Kappa) %*% (L);
             #- mahalanobis(L, 0, mNIW.Kappa);  #<-- more compact
             
    mNIW.Psi <- Psi + S + C;
  }
  
  # finally, do the heavy lifting given these posterior parameters
  result = rMNIW.Rcpp(m, mNIW.Mu, mNIW.Kappa, mNIW.Psi, mNIW.df)
  
  # rename the posterior samples to match the API
  names(result)[names(result) == "X"] = "B" 
  
  return(result)
}


test.rMNIW <- function() { #TODO: move this somewhere more appropriate
  message("Smoketesting MNIW.Rcpp")
  print(rMNIW.Rcpp(n, matrix(1:6,2,3), diag(2), diag(3), 10))
  message("---------------------------------")
}

#test.rMNIW()

