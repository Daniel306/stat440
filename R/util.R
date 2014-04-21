# util.R
# utility functions created in the course of doing the NIW project.
# many of these might already exist in R; please call us on our oversights.
#
# TODO: the docstrings need to be fleshed out
# TODO: these functions need unit tests

is.square <- function(M) {
   # test if a matrix is square
    length(dim(M)) == 2 && (nrow(M) == ncol(M))
}

is.symmetric <- function(M) {
   #test is a matrix is symmetric
    is.square(M) && all(t(M) == M)
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



marginalize <- function(f, arity, dims) {
  # numerically marginalize over a chosen set of inputs to a function
  # 
  # this is somewhat like currying, but it involves a call to integrate for each dim you marginalize over
  # e.g. marginalize(f, 7, c(2,3)) -> f(1,4,5,6,7)
  #   marginalize(f, 7, c(1,3,5,6,7)) -> f(2,4)
  #  marginalize(dnorm, 1, 1) -> f() == 1 (#because PDFs must sum to 1)
  # f needs to be vectorized (because integrate wants it to be vectorized). If you cannot do
  # so f must be able to take a MATRIX of (n, arity)
  #  -> the return value will also be vectorized
  # TODO: clean up this docstring
  # TODO: allow specifying the endpoints for each dim (instead of assuming c(-Inf, +Inf))
  # TODO: this function is extremely slow, and probably also numerically unstable
  #
  # returns:
  #   a function whose arity is arity-length(dims)
  #   (actually, a function which takes a single vector of that length)
  #   (TODO: figure out how to actually handle v-arity functions in R)
  # TODO: support dims being given as negatives (as apply() and indexing)
  
  # TODO: is there some kind of magical optimization which will allow us to precompute the integrate()s? Redoing length(dims) integrates at EVERY call it extremely painful
  
  if(length(dims) == 0) {
    # base case: we're already done marginalizing
    f
  } else { 
    # recursive case: take the first of dims off and marginalize over the rest, then integrate()
    function(v) {
      # vectorization is a jerk
      # it means that v needs to be a matrix...
      v = matrix(v) #coerce the input to a matrix..
      n = dim(v)[1]
      d = dim(v)[2]
      
      stopifnot(d == arity - length(dims)) #arity of the marginalized function must be arity-length(dims)
      
      # curry f over the arguments we know about
      c = function(cv) {
        # Note: no typechecks here because this is an internal function and we are (supposed to be) careful about what we pass it
        #
        # --> cv and v are both matrices
        # cv holds the variates that we are marginalizing over
        # v holds the variates we are pushing up to the caller
        # merging them into a single thingy
        # the way we do this is by cbind()ing and then rearranging the columns using R's indexing magic
        #
      }
      
      # recurse!
      g = marginalize(c, arity, dims[-1])
      
      # g needs to be vectorized for integrate to work on it
      # which means the output of marginalize() needs to be vectorized
      #(which is probably a good idea)
      #(...but that means that v 
      # i need to write 
      # # bear with my weird notation for a moment:
      # say f is f(a,b,c,d). then
      #  marginalize(f, (b,c)) is
      #  marginalize(marginalize(f, b), c) #... is this true?
      # 
      return(integrate(g, -Inf, +Inf))
    }
  }
}

# sketchy test

dmvnorm <- function(X, Mu, V, log=FALSE) {
  if(log) {
    stop("log probabilities are not supported yet")
  }
  #TODO: typechecks
  d = length(X) 
  C = sqrt((2*pi*det(V))^d) #TODO: fix this scaling constant
  exp(-mahalanobis(X, Mu, V)/ 2) / C
}

test.marginalize <- function() {
  # what's an easy trivariate function we can look at?
  #.. how about dmvnorm
  source("test.constants.R")

  # marginalize requires currying over any distributional parameters
  # whatvever, we can use Psi instead of V
  f = function(X) { dmvnorm(X, kMu, kPsi) }
  # perform the tricky part 
  p = marginalize(f, 3, c(2,3))

  # perform the *hard* part
  x = seq(-100, 100, length.out=55)
  plot(x, p(x), main="Marginalization test")
}
test.marginalize()