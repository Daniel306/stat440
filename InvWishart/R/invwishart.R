
# notation: ~ 'distributed as'
#           W 'Wishart distribution'
#       inv(W)'inverse Wishart distribution' -- this is the distribution of matrices whose inverses are Wishart distributed
#           V 'a variance matrix'
#           n 'an int, the 'degrees of freedom' in question" (see [Chi-Squared Distribution]] and [[Wishart Distribution]] for defitinion of the last two)
# Wikipedia says
# X ~ W(V,n)  <=>   inv(X) ~ inv(W)(inv(V), n)
# so, a trivial (slow, but correct) algorithm for (approximating) inverse wisharts is:

require("matrixcalc");

multigamma <- function(x, p) {
  
}

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
  
  pWishart(solve(X), solve(V), df);
}

# multivariate cdfs don't make much sense, because you would
# need to specify what direction you're integrating in
# and more critically because they need to take a ~boundary~ to integrate up to, they can't just take a single ~point~.
# so we don't provide one
# you can use integrate() to approximate a region of interest.

rInvWishart <- function(n, df, Sigma) {
  rWishart(n, df, solve(Sigma));
}

# tests:

