# util.R
# utility functions created in the course of doing the NIW project.
# many of these might already exist in R; please call us on our oversights.
#
# TODO: the docstrings need to be fleshed out
# TODO: these functions need unit tests

source("util.math.R")

# TODO: multidimensionaltranspose
# useful for correcting apply() when it gets it wrong
#  - 1) vectorize the matrix
#  - 2) make clever use of a recursive rep() call to generate the order of the elements (eg. map (k, i, j) -> (i, j, k))
#  - 3) call order() on that order to compute the inverse permutation
#  - 4) shuffle the vectorized matrix by indexing the vectorized matrix by that permutation
#  - 5) put dimensions back on it

..getitem.. <- function(S, idx) {
  # wrapper to work around that in R you cannot say 
  # extracts elements "idx" from structure S
  # an element in idx being NULL is like the syntax  the syntax S[,,3]
  # BUG: this will NOT work if idx contains subvectors like S[1:3, 2:5] (ie idx is a matrix..)
  # TODO: maybe if we can just figure out how to create symbols like R does when you say S[2,,3] then we do can just use that instead of this gunky
  
  #message("..getitem..")
  #print(dim(S))
  #print(idx)
  
  # R magic!: this line courtesy of mrflick in irc.freenode.net/#R
  #do.call(`[[`, c(S, as.list(idx))) #[[ does everything [ does and more
  # ^ but this is finicky if you skip / don't define some dims

  # but there's slice.index which can help:
  # we convert idx into a selection mask
  stopifnot(length(idx) == length(dim(S)))
  mask = rep(TRUE, prod(dim(S)))
  dim(mask) == dim(S)
  for(i in 1:length(dim(S))) { 
    if(is.na(idx[i])) next; # irrelevant; skip it
    mask = mask & slice.index(S, i) == idx[i]
  } #warning! indexing by NA returns NA! this is probably not what you want!!
  stopifnot(all(!is.na(mask)))
  S[mask]
}


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

reasonable_subset <- function(D, length.out=1000) {
  # plotting too many points causes lag
  # reasonable_subset evenly reduces the number of samples in D evenly
  # (in other places, this operation is called decimation(TODO: FACTCHECK))
  #
  # returns:
  #  the dataframe reduced to length.out or its original length, whichever is smaller.

  n = dim(D)[1]
  idx = intseq(1, n, length.out=length.out)
  D[idx,]  
}






marginalize <- function(f, arity, dims, lower=-Inf, upper=+Inf, ...) {
  # numerically marginalize some chosen set of variates out of a pdf
  # 
  # this is somewhat like http://en.wikipedia.org/wiki/Partial_application, but that 
  # each variable is totally removed by an integrate() instead of fixed to a single point.
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
      #print(order(perm)) #DEBUG
      V = V[order(perm)] #the inverse of a permutation is its reverse
      
      # finally, call the original function on the curried data cv mixed with the passed data v
      stopifnot(length(V) == arity)
      #print(V) #DEBUG
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
#test.marginalize()





marginals_do <- function(M, f) { #TODO: pull out into util
  # loop over all the marginals of MultiS M
  # precondition: M's dim is (a,b,c,d,...,n); marginals are taken of that last dimension, n    
  # precondition: f is f(idx, v) where idx will be a vector containing the length(dim)-1 indexes and v is the marginal vector
  # postcondition: nothing; you need to do every with side-effects in f
  # loop over all interesting_dimensions and plot each marginal
  # recursive something someting  
  # This is something like the sum of map() and enumerate() (which are not R functions)
  #   but it is dimensionality aware.
  # TODO: give a flag that marks M as symmetric, so that this only computes on the upper triangle
  #       this is a very-very-special-case flag, though... the cleaner (but inefficient) solution might be to force the lower triangle to NA before running marginals_do
  
  kept_dimensions = 1:(length(dim(M))-1) #drop the last dimension
    
  R<-function(idx, dims) {
  
    if(length(dims)==0) {
    
      # base case
      # we've bottomed out and have an actual index in hand
      # so use it
      #message("bottomed out at idx=")
      #print(idx) #DEBUG

      idx = c(idx,NA) #this has to be extended by one NA in order to extract the vector we care about; see ..getitem..
      v = ..getitem..(M, idx)
      #message("WE EXTRACTED THIS") #DEBUG
      #print(v)
      
      v = as.vector(v) # just in case (i think if idx == c(), which should correspond to M[], we will *not* get a vector out: R will preserve the dimensions (but only in that case (which shouldn't be possible, if the recursion is correct))
      f(idx,  v)
      
    } else {
    
      # recursive case
      # split the top dimensions from the body
      # note how we extract the *number* of elements in the dimension d from the index of the dimension itself dims[1]
      d = dim(M)[dims[1]]; dims = dims[-1];
      for(i in 1:d) {
        # recurse!
        R(append(idx, i), dims)
      }
      
    }
    
  }
  R(c(), kept_dimensions) #kick off the recursive loop  
  return(NULL) #don't return whatever R() happens to construct; that would be terrible
}


idx2str <- function(idx) {
  # convert an index vector as passed from marginals_do
  # into a string
  # c()      -> ""
  # c(1,2,3) -> "[1,2,3"]
  # c(1)     -> "[1]"
  idx[is.na(idx)] = ""
  idx = paste(idx, collapse=",") #could also do this with a complicated do.call, but paste() covers this case helpfully for us with 'collapse'
  if(idx == "") { #special case: no indices means all the indecies #XXX untested
    return(idx)
  } else {
    return(paste("[", idx, "]", sep="")) #but otherwise wrap the indecies in square brackets
  }
}

