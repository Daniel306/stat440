# test.MultivariableRegression.R
#
#
source("MultivariableRegression.R")


# originally n = 27, m = 1e5
test.EversonMorris <- function(n=27, m=1e7) {
  # Smoketest for lm.multivariable
  #
  # simulate multivariable normal test data
  # then see if the bayesian fitter can pull out the
  # correct (artificial and frequentist) coefficients.
  # Parameters based on Everson & Morris [2000]
  #
  # Note: the proper Everson Morris paper used a *random effects* model which
  #  turns out(?) to be algebraically equivalent to a multivariable model.
  #  This function only uses the multivariable model directly.
  #
  # args:
  #  n: number of observed samples to take
  #     the default is purposely small, to reflect
  #     a realistic data-gathering situation.
  #  m: number of posterior samples to take
  #
  # returns:
  #  nothing; instead, results are printed as work is done.

  d = 1;
  q = 2;
  A = matrix(c(3.38,-.77,-.77,2.55), q, q)
  B = matrix(c(0,0), d, q)
	
  message("True coefficients")
  print(B)
  message("True covariance")
  print(A)
  
  data = rmultivariableregression(n, B, A);
  
  message("Hiding true values from ourselves")
  rm(B, A, d, q)
  
  message("Data is:")
  message("Y:")
  print(data$Y)
  message("X:")
  print(data$X)

  X.sq <- crossprod(data$X);
  beta.hat1 <- solve(X.sq, t(data$X) %*% data$Y);
  message("beta hat point estimate")
  print(beta.hat1)
  
  message("Sampling posterior distribution")
  #  m = 2 #DEBUG
  result <- lm.multivariable(m, data$X, data$Y)
 
  d = dim(result$B)[1]
  q = dim(result$V)[1]
  
  message("Recovered dimensionality: q = ", q, " and d = ", d)
  message("Estimated B")
  B.hat = apply(result$B, -3, mean) #take mean()s across the 3rd dimension;
                            # results in a q x d matrix (after unflattening apply's results)
  dim(B.hat) = c(q,d)
  print(B.hat)
  message() #newline

  message("Estimated V")
  V.hat = apply(result$V, -3, mean) #take mean()s across the 3rd dimension;
                            # results in a q x q matrix (or it would, but apply flattens its results)
  dim(V.hat) = c(q,q)
  print(V.hat)
  message() #newline
  
  # TODO: compute (with 'quantile') the confidence intervals for each B and V
  #
  # additionally, we have a whole sample from which we can do bootstrap-like things, compute functions of the data, etc
  # but for this simple test, getting the right coefficients is good enough.
  list(V = result$V, B = result$B, X = data$X, Y = data$Y)
}

