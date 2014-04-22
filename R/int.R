
vectorize <- function(f) {
  # vectorize a univariate function
  #
  # this uses sapply(), so it's reliable and as slow as doing it by hand
  # it can make your code more readable, though.
  
  function(x) {
    sapply(x, function(x) { f(x) })
  }
}

# multiintegrate
integrate.multi <- function(f, lower, upper, ...) {
  # perform a k-ary, hyperrectangular, integral, numerically
  # this is not API compatible with integrate() because it only returns a single number
  # if any of the subintegrations fails the whole thing will probably just blow up
  #
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

# test
source("util.R")
dmvnorm <- function(X, Mu, V, log=FALSE) {
  if(log) {
    stop("log probabilities are not supported yet")
  }
  
  # typechecks
  stopifnot(is.vector(X) && is.vector(Mu) && is.matrix(V))
  d = length(X) 
  stopifnot(length(Mu) == d)
  stopifnot(is.symmetric(V))
  stopifnot(dim(V)[1] == d)
  
  C = sqrt((2*pi)^d*det(V))
  exp(-mahalanobis(X, Mu, V)/ 2) / C
}
# ^warning: dmvnorm is *not* vectorized

source("test.constants.R")
r = integrate.multi(function(X) {
  dmvnorm(X, kMu[1:2], kPsi[1:2, 1:2])
}, c(-Inf, -Inf), c(Inf, +Inf))
message("integral came out to: ", r)