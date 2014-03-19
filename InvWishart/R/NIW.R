
require("matrixcalc")
require(MASS); #for rmvnorm
rmvnorm <- mvrnorm #why
mvrnorm <- NULL


rInvWishart <- function(n, df, Sigma) {
 # inverse wishart sampler
 message("rInvWishart")
 message("----------------------")

 print(Sigma)
 print(class(Sigma))
 print(is.square.matrix(Sigma))
 
 W = rWishart(n, df, solve(Sigma));
 for(i in 1:n) {
   W[,,i] = solve(W[,,i]) 
 }
 W
}

rNIW.naivesample <- function(Mu, kappa, Psi, df) {
  # naive slow-as-mud Normal-Inverse-Wishart sampler
  print(Psi)
  V = rInvWishart(1, df, Psi)[,,1];
  message("V is")
  print(V)
  X = rmvnorm(1, Mu, V/kappa);
  print(X)

  list(X=X, V=V);
}

rNIW.snappysample <- function(Mu, kappa, Psi, df) {

# attempt to write out the efficient algorithm (but in R, of course)
     # construct L_i
     # construct W_i (?)
     # construct U_i (?)
     # construct V_i

list(X=X, V=V);
}

#(one thing this setup doesn't do is precompute inv(Psi) ahead of time)

# 

rNIW.sample = rNIW.naivesample

rNIW.wrapper <- function(n, Mu, kappa, Psi, df) {
  message("hurp")
  # n: the desired number of sample points (X, V) ~ NIW(Mu, kappa, Psi, df)
  # df need not be integer; instead it can be anything in the domain of the gamma function
  
  # this subroutine is just the safe(ish) type-checked outer call
  # which calls rNIW.sample to get the real samples (so you can experiment with swapping out different samplers)
  
  Mu = as.vector(Mu)
  Psi = as.matrix(Psi)
  
  d = dim(Psi)[1]
  stopifnot(d == dim(Psi)[0]) #Psi must be square
  stopifnot(d == length(as.vector(Mu)))  #Mu must be the same dimensionality as the covariance matrix it supposedly goes with
  
  # pre-allocate space
  X = rep(0, d*n);
  dim(X) = c(d, n);
  
  V = rep(0, d^2*n);
  dim(V) = c(d,d,n);
  
  # generate the n samples
  for(i in 1:n) { 
     message(i)
     sample = rNIW.sample(Mu, kappa, Psi, df);
     X[,i] = sample$X
     V[,,i] = sample$V
  }
  list(V=V, X=X)
}


rNIW = rNIW.wrapper

kMu = c(0.25, -1.5, 0.33)
kKappa = 1
kV  = cbind(c(2.84, 0.43, 0.16), c(0.43, 1.52, -0.24), c(.16, -.24, 4.49))
kDF = 3.3


Z = rNIW(6, kMu, kKappa, kV, kDF)

print(Z)
