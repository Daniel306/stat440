
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
integrate.multi(f, left, right, ...) {
  # perform a k-ary, hyperrectangular, integral, numerically
  #
  # TODO: clean up the 'extra args' part
  # ...: extra args to integrate()
  # note: this i

  # Typechecks
  stopifnot(is.vector(left) && is.vector(right)) #note: scalars are vectors, in R (but vectors are not matrices!)
  p = length(right)
  stopifnot(p == length(left))
  stopifnot(p > 0) #TODO: is it meaningful to integrate over a function with no arguments?
  #          -> ..sort of. the result is a rectangle: f() * (right-left)
  #  except that isn't right, because, as written, a function with no args is marked by right and left being empty
  #  so what does that mean?
  # that there are no
  # still, it seems like this can be written shorter if we define a p==0 case..
  
  if(p == 1) {
    # base case
    integrate(f, left, right, ...)
  } else {
    # recursive case
    # note that the function we pass down must be able to handle vector args
    integrate(function(x) {
      # integrate needs to take a vectorized function
      # doing that with matrices is tricky
      # so for now, do it with sapply(), as inefficient as that is
      sapply(x, function(x) {
        # in here, x is scalar
        # recurse!
        integrate.multi(function(v) {
          # do something...
          #???
          # Splay scalar x against all the values in v
          # we need a function which computes
        }
        , left=left[-1], right=right[-1])
      })
    }, left=left[1], right=right[1], ...)
  }
}

# test
dmvnorm <- function(X, Mu, V, log=FALSE) {
  if(log) {
    stop("log probabilities are not supported yet")
  }
  #TODO: typechecks
  d = length(X) 
  C = sqrt((2*pi*det(V))^d) #TODO: fix this scaling constant
  exp(-mahalanobis(X, Mu, V)/ 2) / C
}
# ^warning: dmvnorm is *not* vectorized

source("test.constants.R")
r = intergrate.multi(function(X) {
  dmvnorm(X, kMu, kPsi)
}, c(-Inf, -Inf, -Inf), c(0, 0, +Inf))
message("integral came out to: ", r)