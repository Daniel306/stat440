# util.R
# utility functions created in the course of doing the NIW project.
# many of these might already exist in R; please call us on our oversights.
#
# TODO: the docstrings need to be fleshed out
# TODO: these functions need unit tests

is.square <- function(M) {
    # test if a matrix is square
    # TODO: this appears to exist as matrixcalc::is.square.matrix already??
    length(dim(M)) == 2 && (nrow(M) == ncol(M))
}

is.symmetric <- function(M) {
    #test if a matrix is symmetric
    is.square(M) && max(t(M) - M) < 1e-6
}


vectorize <- function(f) {
  # vectorize a univariate function
  #
  # this uses sapply(), so it's reliable and as slow as doing it by hand
  # it can make your code more readable, though.
  # TODO: allow x to be a matrix
  # TODO: compare to Vectorize() and maybe toss this one
  
  function(x) {
    sapply(x, function(x) { f(x) })
  }
}



clip <- function(x, range) {
  # hard-clip values to a specified range
  #
  # args:
  #  x:     the value to clip
  #  range: a pair c(l, u), as in the output of the range() function.
  #         Use l = -Inf or u = +Inf to only partially clip.
  
  # TODO: typechecks
  l = range[1]
  u = range[2]
  
  max(l, min(u, x))
}

intseq <- function(from=1, to=NULL, by=NULL, length.out=NULL) {
  # construct a list of indecies
  # 
  # this is like seq() except that 'seq' doesn't guarantee its results can be used as indecies.
  # this is NOT the same as seq.int, which seems mostly to just be an alias for seq
  # usage: intseq(n) -> 1:n
  #        intseq(m,n) -> m:n
  #        intseq(m,n,by=k) -> c(m, m+k, m+2k, ..., n-ish) where 'ish' is the remainder of (n-m)/k
  #
  # args:
  #  TODO
  #
  # returns:
  #   a vector of ordered, positive integers
  #
  # TODO: test corner cases; as written this might *sometimes* have off-by-one bugs where the output length doesn't match length.out

  if(missing(to)) {
    to = from;
    from = 1;
  }

  if(missing(by)) {
    if(missing(length.out)) {
      by = 1;
    }
    if(length.out > (to-from)) {
      by = 1;
    } else {
      by = floor((to - from) / (length.out - 1)) #this formula taken from help(seq); I am unsure why the -1 is necessary, but it is demonstrably.
    }
  }

    
  if(!(by==round(by))) {
    stop("assertion: 'by' should be integer, but by=",by," and round(by)=", round(by))
  }
  stopifnot(by > 0)
  
  seq(from, to, by)
}


cummean <- function(v) {
  # compute the "cumulative mean" of a vector: considering each entry as a new sample point, give the sample mean up to that point
  # 
  # args:
  #  v: a numeric vector
  #
  # returns:
  #   a vector the same length as v
  
  cumsum(v) / (1:length(v))
}
 
dmvnorm <- function(X, Mu, V, log=FALSE) {
  # the multivariate normal density function (pdf)
  # warning: dmvnorm is *not* vectorized
  # TODO: this function isn't really a utility; it belongs ..somewhere else
  # TODO: docstring
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


marginalize <- function(f, arity, dims, lower=-Inf, upper=+Inf, ...) {
  # numerically marginalize some chosen set of variates out of a pdf
  # 
  # this is somewhat like currying, but each curried variable corresponds to an integrate()
  #
  # the marginalized function is EXTREMELY slow becaues a every evaluation of it
  #  requires a length(dim)-deep stack  with at least 10, probably more  (however many integrate() decides) 
  #  forks at each branch.
  # It is also prone to numerical instability.
  # Basically: the curse of dimensionality bites this function and there's no (generic, numerical) way around that.
  #
  # args:
  #   f: the function to marginalize over; f should take a single argument, a vector, and return a single scalar
  #   arity: the number of variates that f takes; that is, the length of the single vector argument to f.
  #         be careful to specify this correctly. marginalize() has no way of typechecking that you are giving the correct
  #         or that f is interpreting the vector marginalize ends up handing it correctly
  #   dims: a vector containing which variates to marginalize out
  #  lower, upper: vectors giving the boundaries of dims to marginalize over
  #    this can speed up evaluation, but possibly at the expense of correctness.
  #    You must know something about your distribution to use these correctly.
  #    These should be the same length as the number of dimensions.
  #       as a special case, if either is scalar, it is used as the boundary for ALL dimensions.
  #  ...: extra arguments to 'integrate()'
  #
  # e.g. if you have a function f(a,b,c,d,e,f,g) then marginalize(f, 7, c(2,3)) -> f(a,d,e)
  #   marginalize(f, 7, c(1,3,5,6,7)) -> f(2,4)
  #  marginalize(dnorm, 1, 1) -> f() == 1 (#because PDFs must sum to 1)
  #
  # returns:
  #   a function which takes a single vector of ariity
  #
  # TODO: figure out how to actually handle v-arity functions nicely in R
  # TODO: support dims being given as negatives (as apply() and indexing)
  # TODO: is there some kind of magical probability theory trick which allows precomputing the integrate()s?
  #             Redoing length(dims) integrates at EVERY call it extremely painful
  #             ESPECIALLY because integrate.multi can't handle vectorized calls
   
   # Typechecks
  stopifnot(unique(dims) == dims)
  stopifnot(all(1 <= dims && dims <= arity))
  # By the pigeonhole principle, these two lines also imply this:
  #stopifnot(length(dims) <= arity)

  function(v) {
    #Typechecks
    stopifnot(is.vector(v))
    p = length(v)
    q = length(dims)
      
    stopifnot(p + q == arity) #arity of the marginalized function must be the number of marginalized variables (length(dims)) and the number of variates given to us here
    # above we checked that length(dims) <= arity,   so   arity - length(dims) >= 0 
      
      
    # if the integration bounds are scalar,
    # vectorize them
    if(length(lower) == 1) {
      lower = rep(lower, q)
    }
    if(length(upper) == 1) {
      upper = rep(upper, q)
    }
    stopifnot(length(lower) == length(upper))
    stopifnot(length(lower) == q)
    
      
    # here's the idea: to marginalize, first curry the function so that it is only a function of the marginalized variates
    # , then use integrate.multi to eat up the marginalized variates
    cr = function(cv) {
        
      # Note: no typechecks here because this is an internal function and
      # we are (supposed to be) careful about what we pass it
       
      # construct a single
      V = c(cv, v)
      # but V is in the wrong order: the curried arguments cv need to go where 'dims' says they should
      # and the others v need to go in the remaining spots
      # We can think of this problem as that V has been permuted from the true order, and we want to undo the permutation
      # In (pseudo-)R notation, our permutation is
      #   c(dims, -dims)
      # For example:
      #   c(2,4,1,3) means 
      #  column 2 --> 1
      #  column 4 --> 2
      #  column 1 --> 3
      #  column 3 --> 4
      # but we want the inverse, which would be c(3,1,4,2)
      # it turns out that for some cases, rev() is enough
      # but in general, we actually can (mis)use order() to 
      # TODO: factor this difficult code (and the reasoning behind it) to a subroutine
      cv_dims = dims
      v_dims = (1:arity)[-dims] #R doesn't let us directly mix negative indexing with positive, so we need to do this
      perm = c(cv_dims, v_dims)
      print(order(perm))
      V = V[order(perm)] #the inverse of a permutation is its reverse
      
      # finally, call the original function on the curried data cv mixed with the passed data v
      stopifnot(length(V) == arity)
      print(V) #DEBUG
      f(V)
    }
    
    # now that we have curried, marginalize over everything else
    integrate.multi(cr, lower=lower, upper=upper, ...)
  }
}

# sketchy test
test.marginalize <- function() {
  # what's an easy trivariate function we can look at?
  #.. how about dmvnorm!
  # warning: this function is extremely slow. Call with time to spare.
  source("test.constants.R")
  kV = solve(kPsi)
  #kMu = kMu[1:3] #the constants are 4 dimensional
  #kV = solve(kPsi[1:3,1:3]) #but that's too much to test.

  d = length(kMu)  
  # marginalize requires currying over any distributional parameters
  # whatvever, we can use Psi instead of V
  # the marginal of a multinormal should be normal
  f = function(X) { dmvnorm(X, kMu, kV) }
  # perform the tricky part 

  margin = 2
  p = marginalize(f, d, (1:d)[-margin], rel.tol=1e-3)
  
  # perform the *hard* part
  x = seq(-3, 3, length.out=9)
  plot(x, vectorize(p)(x), main="Marginalization test")
  
  # the expected marginal is just found by simply dropping
  # the irrelevant dimensions from the parameters
  message("Multivariate normal distribution:")
  message("mu:")
  print(kMu)
  message("V:")
  print(kV)
  
  for(margin in margin) {
    abline(v=kMu[margin], lty="dashed", col=margin)
    dist = function(x) { dnorm(x, kMu[margin], kV[margin,margin]) }
    lines(x, dist(x), lty="dashed", col=margin) #plot the expected density
  }
}
test.marginalize()
