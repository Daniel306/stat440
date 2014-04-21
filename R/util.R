# util.R
# utility functions created in the course
# many of these might already exist in R; please call us on our oversights
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
    } else {
      by = (to - from) / (length.out - 1) #this formula taken from help(seq); I am unsure why the -1 is necessary, but it is demonstrably.
    }
  }
    
  stopifnot(by==round(by)) #for 
  
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

