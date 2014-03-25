

# Style plan:
#  -follow Google's R style guide

# API plan:
#  -do we conform to the API of rWishart ("(df, V)") or use the dist'n parameters ("(V, df)")?


# Testing plan:
#  -testing is ...hard?
#  - well, for one, we can generate rWisharts, then generate their Inverses
#  - use integrate() to check summing to 1 (tho as d increases you will)
# correctness testing:
#  - we need to compare distributions; this is hard, because our distributions are on matrices (or generally, vectors)
#  - we can snitch HW1's code to make plots showing marginal histograms against kernel density estimates for
#     a vectorized (might need some retrofitting to generalize it to all vectors? or maybe we just insist its for matrices)
#  - 
#
# performance testing
#  - compare the naive R implementation with an Rcpp implementation
#  - compare the naive SciPy implementation with an
#  - compare the less naive implementation (ie the one given us by Lysy) in: pure R, Rcpp, Scipy, Pyrex
#  -
#  comparison algorithm is:
#   - parameter sweep over (n, d, df, implementation) -- use a random V in each case, and record runtimes
#   - plot runtimes as functions of each dimension
#   - plot runtimes in 5d space via magic (actually, break it down first splitting by implementation, then
#   - include confidence intervals (bands??)


# notation: ~ 'distributed as'
#           W 'Wishart distribution'
#       inv(W)'inverse Wishart distribution' -- this is the distribution of matrices whose inverses are Wishart distributed
#           V 'a variance matrix'
#           n 'an int, the 'degrees of freedom' in question" (see [Chi-Squared Distribution]] and [[Wishart Distribution]] for defitinion of the last two)
#

#######################################3
### IMPORTS
###

require("matrixcalc");

#######################################
### UTILS

multigamma <- function(x, d) {
  # usage: multigamma(x, d)
  #
  # x: the point to evaluate at
  # d: the dimensionality of the multigamma
  # 
  prod(x + (1-1:d)/2) / sqrt(sqrt(pi^(d*(d-1))))
}


multigamma <- function(d) {
  # usage: multigamma(d)(x)
  
  function(x) {
  #
  # x: the point to evaluate at
  # d: the dimensionality of the multigamma
  # 
  prod(x + (1-1:d)/2) / sqrt(sqrt(pi^(d*(d-1))))
  }
}


#######################################
### PDFS

pWishart <- function(X, V, df) {
  #TODO: docs
  #TODO: tests
  
  if (!is.square.matrix(X))  {
    stop("Wishart matrices must be square");
  }
  
  if (!is.square.matrix(V))  {
    stop("Covariance matrices must be square");
  }
  
  d = dim(X)[0]  #dimensionality of the matrix: note the matrix is actually dxd, but it comes from an associated d-dimensional Normal (sampled df times)
  if(!(d == dim(V)[1] && d == dim(V)[2])) { #XXX this is redundant: is.square.matrix already checked..
    stop("Wishart matrices must have the same dimension as their generating covariance matrix");
  }
  
  if(df < d) {
    stop("Wishart matrices cannot have dimensionality less than their degrees of freedom");
  }
  
  #TODO:
  # enforce positive-definiteness on V
  # enforce matrix-ness
  # TODO: does X also have to be positive definite? (no.. it's random..?)
  
  # from wikipedia:
  # rewritten a bit to reduce the sqrt to once
  # the matrix inversion (solve(V)) is dangerous
  # XXX can trace(inv(V)*X) be reduced?
  # ..well, trace(inv(V)*X) = trace(inv(L'L)*X) = trace(inv(L)*inv(L')*X) = {rotation of tr, maybe slightly wrong} trace(inv(L')*X*inv(L)) = ???
  sqrt(exp(-matrix.trace(solve(V)%*%X) ) * det(X)^(df - d - 1) / det(V)^df / 2^(df*d))  / multigamma(p/2,d)
}

pInvWishart <- function(X, V, df) {
  # inverse wishart PDF
  # 
  # args:
  #   X: a 
  #   V: a variance-covariance matrix
  #   Sigma: the number of 


  # so, say we know how to compute P(W(V,n) = X)
  # from that, we can compute P(inv(W)(Q, n) = M) because if M ~ inv(W)(Q, n) then inv(M) ~ W(inv(Q), n) so we ask pWishart(inv(M), inv(Q), n)
  # TODO: derive this more clearly/formally/in latex-ly
  #  ---like, perhaps it is true that these two dists are equidistributed, but does that mean that the PDFs are the same??
  
  pWishart(solve(X), solve(V), df); #XXX FIXME: this needs to be scaled by the jacobian of V^-1. See Sawyer-Wishart.pdf.
}


#######################################
### CDFS

# multivariate cdfs don't make much sense, because you would
# need to specify what direction you're integrating in
# and more critically because they need to take a ~boundary~ to integrate up to, they can't just take a single ~point~.
# so we don't provide one
# you can use integrate() to approximate a region of interest.


#######################################
### SIMULATORS

rInvWishart.naive <- function(n, df, Sigma) {
  # naive algorithm:
  # Wikipedia says
  # X ~ W(V,n)  <=>   inv(X) ~ inv(W)(inv(V), n)
  # so, a trivial and correct (but slow?? TODO: TEST) algorithm for (approximating) inverse wisharts is:
  
  rWishart(n, df, solve(Sigma));
}


rInvWishart.extremelynaive <- function(n, df, Sigma) {
  # extremely naive algorithm:
  # the Inverse Wishart distribution is defined as the distribution of matrices who are the inverse of matrices that are Wishart distributed
  
  R = rWishart(n, df, Sigma);
  # now, invert all of them
  for(i in 1:n) {
     R[,,i] = solve(R[,,i])
  }
  R
}


rInvWishart = rInvWishart.naive
