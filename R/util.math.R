# util.math.R
#
# heavy-lifting math


dmvnorm <- function(X, Mu, V, log=FALSE) {
  # the multivariate normal density function (pdf)
  # X can be a single vector, or a matrix (n, d) (<-- note, this contradicts the convention used elsewhere in this project)
  # TODO: this function isn't really a utility; it belongs ..somewhere else
  # TODO: docstring
  # matching rmvnorm in rNIW.cpp
  
  # typechecks
  stopifnot(is.vector(Mu) && is.symmetric(V))
  d = length(X) 
  stopifnot(length(Mu) == d)
  stopifnot(is.symmetric(V))
  stopifnot(dim(V)[1] == d)
  if(is.vector(X)) {
    X = t(matrix(X)) #a vector is a single sample and mahalanobis expects each sample to run across columns a row
  }
  stopifnot(dim(X)[2] == d)

  # in normal space (this formula from Wikipedia)
  #C = sqrt((2*pi)^d*det(V))
  #exp(-mahalanobis(X, Mu, V)/ 2) / C
  
  # in logspace
  c = (d*(log(2) + log(pi)) + log(det(V))) #normalizing constant
  p = -mahalanobis(X, Mu, V) #<-- this vectorizes. X needs to be (n, d)
  p = (p - c)/2;
  
  if(!log) {
    p = exp(p)
  }

  return(p)
}

trace <- function(M) {
  # trace of a matrix
  # TODO: typechecks
  # Warning: This shadows the trace function meant for debugging
  sum(diag(M))
}

# multiintegrate
integrate.multi <- function(f, lower, upper, ...) {
  # perform a k-ary, hyperrectangular, integral, numerically
  # this is not API compatible with integrate() because it only returns a single number
  # if any of the subintegrations fails the whole thing will probably just blow up
  ## f should be a function taking a single vector argument, which is all the variates that it is a function of
  # 
  # TODO: support vectorized fs, which take a whole set of points (ie an n x p matrix) to sample at
  # TODO: clean up the 'extra args' part
  # ...: extra args to integrate()
  # note: this i

  # Typechecks
  stopifnot(is.vector(lower) && is.vector(upper)) #note: scalars are vectors, in R (but vectors are not matrices!)
  p = length(upper)
  stopifnot(p == length(lower))
  stopifnot(p > 0) #TODO: is it meaningful to integrate over a function with no arguments?
  #          -> ..sort of. the result is a rectangle: f() * (right-left)
  #  except that isn't right, because, as written, a function with no args is marked by right and left being empty
  #  so what does that mean?
  # that there are no
  # still, it seems like this can be written shorter if we define a p==0 case..
  
  g = 
  if(p == 1) {
    #base case
    f
  } else {           
    #recursive case
    function(x) {   #<-- this x is a scalar, the first argument to f
      o = function(v) { #<-- this v is a vector *of* scalars: the rest of the arguments
                     # though v is a vector, neither of these functions
                     # are vectorized (in the sense of mapping multiple
                     # input points to an equal number of output points)
                     # (that happens later)
        f(c(x, v))
      }
    integrate.multi(o, lower=lower[-1], upper=upper[-1]) #-1 deletes the 1st element from the vectors
    }
  }
  
  integrate(vectorize(g), lower=lower[1], upper=upper[1], ...)$value
}


test.integrate.multi <- function() {
  source("test.constants.R")
  r = integrate.multi(function(X) {
    dmvnorm(X, kMu[1:2], kPsi[1:2, 1:2])
  }, c(-Inf, -Inf), c(Inf, +Inf))
  message("integral came out to: ", r, "; (it should be very nearly 1)")
}
#test.integrate.multi()



dinvgamma <- function(x, shape, rate=1, scale=1/rate, log=FALSE) {
  # pdf of the inverse gamma function

  # from wikipedia, of course
  # computation initially done on a logscale
  # lgamma is a different, presumably more accurate, implementation of "ln(gamma(.))"
  
  p = shape*log(scale) - lgamma(shape) - (shape+1)*log(x) - scale/x
  p[x <= 0] = -Inf     # correctly handle x which is out of the support
                       # it would be cleaner to not compute them at all in the first place
                       # but that would be more verbose.
  
  if(!log) {
    p = exp(p);
  }
  
  p
}


# dinvwishart <- function(x,
# already implemented as dIW() in dNIW.R
# ...?
# hmmm