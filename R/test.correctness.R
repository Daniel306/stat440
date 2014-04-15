
#typedef: a multidimensional sample (MultiS) is
#  a matrix M with dim(M) = c(n, c(a,b, c, ...)) 
#  n is the number of samples
#    (a, b, c, ...) is the dimensionality of the object

plot.compare <- (ground, sample, ...) { #XXX NAME
   # plot histograms of the marginals in sample, overlaid with "ground truth" probability density.
   # 
   # ground: the "ground truth"
   #  this can either be a vector of samples, which will be kernel-density estimated, or a pdf function (which must be vectorized!!)
   #  the former is quicker to use: just generate samples in a naive and easy to validate way
   #  But the latter cannot be tainted by mistakes in your naive generator.
   # sample: the "sample"
   #  this is is the sample
   # ...: extra arguments to plot()
   
   if(!is.call(ground)) {
     
     # this approach is a bit convoluted! it is perfectly possible to call lines(density(...)) directly
     # however, doing this special case just feels more clever which is nice
   }

   # display the histogram
   hist(sample, probability=TRUE, ...)
   
   # overlay the pdf
   if(is.call(ground)) {
     lines(sample, ground(sample))
   } else {
     # if caller has not given us a pdf directly, then
     # assume ground is a vector for us to kernel-density estimate
     stopifnot(is.vector(ground) && is.numeric(ground)) #typecheck
     
     lines(density(ground))
   }

  # omitted code, which could be used instead of the special casing in lines() above
    ## replace ground with a (callable!) pdf function 
    ## via kernel density estimation
    # if(!is.call(ground)) {
    # p = density(ground)
    # ground = function(x) {
    #   approx(p$x, p$y, x, rule=2)$y # rule=2 means 'support extrapolation'
    # }
    #}
}


plot.compare.marginals <- function(ground, sample) { #XXX name
  # plot comparison histograms over all the marginals in
  # ground and sample must share the same dimensionality -- all but their first dimension, which is n, the number of samples
  #  and in particular, ground must be a sample itself -- you cannot use pdf functions here

  if(any(dim(ground)[-1] != dim(sample)[-1])) {
    stop("ground and sample have inconsistent dimensions: ", dim(ground)[-1],"  vs ", dim(sample)[-1]) #TODO: this won't print correctly
  }

  # method 1: reshape ground and sample into 2d matrices: n x (number of marginals)
  #   con: you lose the information about ~which~ marginal you're looking at, unless you reconstruct it manually
  # method 2: somehow recurisvely loop down the dimensions, loop over all values in that dimension as you go
  #   con: since [] is a variable arity function call in R, constructing the correct indexing calls will require some R magic
  #        this is probable doable though

  # however
  # e.g. the [1,3] marginal is sample[, 1, 3]: note how every dimension was specified except for the first
  # this doesn't work as written, but it shows the idea of what I want to accomplish: 
  #  plot.compare(ground[, m], sample[, m], main=paste("Marginal ", m))
}

marginals1.do <- function(M, op) {
    # extracts the one-dimensional marginal distributions of MultiS
    #  and performs op(m, v) on them, with m being the index of the marginal (a vector, e.g. c(2,2) means the 2,2th marginal) and v being the vector of values
    #
    # args:
    #   M, a MultiS
    #
    # returns:
    #   a list with entries named "[a,b,c]" -- each
    #    the list is ordered in

    # this is really just a reshape..
        n = dim(M)[1]
            dim(M) = c(n, length(M)/n) #NB: R's quirk: 'length(M)' is the total number of elements in it
                               # length(M) is gauranteed to be divisible by n
                                  }
                                  - 


# 

marginals2.do <- function(M) {

}

#NIW.moment.first
# --> mean?
#NIW.moment.second
# --> ??