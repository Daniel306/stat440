
source("NIW.R")

#typedef: a multidimensional sample (MultiS) is
#  a matrix M with dim(M) = c(n, c(a,b, c, ...)) 
#  n is the number of samples
#    (a, b, c, ...) is the dimensionality of the object

# this is complicated, however, because what is marginalized over does not necessarily come in a sensible matrix shape:
  # for NIW our samples are tuples (X, V) with dim(X)=c(d,1), dim(V) = c(d,d)
  #    we could cram them together so that each sample is a d x (d+1) matrix
  #    but that's awkward and loses pertinent information
  # in general the entries of the tuples may be any size at all


plot.compare <- function(ground, sample, ...) { #XXX NAME
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

plot.NIW.marginals <- function(ground, sample, alg=NULL) {
    # compare (using plot.compare) the marginals of the Normal Inverse Wishart against their expected pdfs
    #
    # args:
    #   ground, sample: lists containing two elements, $X and $V. The first is a c(n, d) matrix, the second is a c(n, d, d) matrix
    #                the each of the n pairs (X[k,], V[k,]) corresponds to one of the Normal Inverse Wishart samples.
    #         ground is assumed to be samples from the true distribution, and creates kernel density estimates
    #         sample has histograms plotted from it
    
    #typechecks
    stopifnot(class(ground) == "list" && class(sample) == "list")
    stopifnot(names(ground) == names(sample))
    for(N in names(ground)) {
      last = length(dim(ground[[N]])) #the length of the last dimension is "n", the number of samples,
                                      # which is irrelevant; we only care that the dimensionality of each sample in ground matches the dimensionality of each sample in sample
      stopifnot(length(dim(ground[[N]])) == length(dim(sample[[N]])))
      if(any(dim(ground[[N]])[-last] != dim(sample[[N]])[-last])) {
        stop("ground$", N, " and sample$", N, " have inconsistent dimensions: ", dim(ground[[N]])[-last],"  vs ", dim(sample[[N]])[-last]) #TODO: this won't print correctly
      }
    }
    
    
    # TODO: factor this into a generic function that takes a list() of matrixables and figures out the labels (eg "V[3,43]") to use automatically 
    # method 1: reshape ground and sample into 2d matrices: n x (number of marginals)
    #   pro: simple
    #   con: you lose the information about ~which~ marginal you're looking at, unless you reconstruct it manually
    # method 2: somehow recurisvely loop down the dimensions, loop over all values in that dimension as you go
    #   con: since [] is a variable arity function call in R (it doesn't just take a tuple like in Python),
    #          constructing the correct indexing calls will require some R call() magic
    #          such magic is probably doable though
  
    # since we know sample came from rNIW, we have stronger conditions than the above:
    # 
    d = dim(ground$X)[1]
    
    # X marginals
    par(mfrow=c(2,2))
    for(i in 1:d) { #the 1st dimension is the 
        plot.compare(ground$X[i,], sample$X[i,], main=paste("X[",i,"]"), sub=alg)
    }
    
    # V marginals
    par(mfrow=c(2,2))
    for(i in 1:d) {
    for(j in 1:d) { #purposely not indented
        plot.compare(ground$V[i,j,], sample$V[i,j,], main=paste("V[",i,",",j,"]"), sub=alg)
    }
    }
}

ks.test.NIW.marginals <- function(ground, sample, alg=NULL) {
    # compare (using plot.compare) the marginals of the Normal Inverse Wishart against their expected pdfs
    #
    # args:
    #   ground, sample: lists containing two elements, $X and $V. The first is a c(n, d) matrix, the second is a c(n, d, d) matrix
    #                the each of the n pairs (X[k,], V[k,]) corresponds to one of the Normal Inverse Wishart samples.
    #         ground is assumed to be samples from the true distribution
    #         sample has histograms plotted from it
    #      for each marginal, ks.test is performed and the results printed
    
    #typechecks
    stopifnot(class(ground) == "list" && class(sample) == "list")
    stopifnot(names(ground) == names(sample))
    for(N in names(ground)) {
      last = length(dim(ground[[N]])) #the length of the last dimension is "n", the number of samples,
                                      # which is irrelevant; we only care that the dimensionality of each sample in ground matches the dimensionality of each sample in sample
      stopifnot(length(dim(ground[[N]])) == length(dim(sample[[N]])))
      if(any(dim(ground[[N]])[-last] != dim(sample[[N]])[-last])) {
        stop("ground$", N, " and sample$", N, " have inconsistent dimensions: ", dim(ground[[N]])[-last],"  vs ", dim(sample[[N]])[-last]) #TODO: this won't print correctly
      }
    }
    # since we know sample came from rNIW, we have stronger conditions than the above:
    # 
    d = dim(ground$X)[1]
    
    # X marginals
    par(mfrow=c(2,2))
    for(i in 1:d) { #the 1st dimension is the 
        p = ks.test(sample$X[i,], ground$X[i,])$p.value #NB: order that ks.test takes its args is swapped from our order
        message(alg," X[",i,"]", ": ", if(p > 0.05) { "same" } else { "different" })
    }
    
    # V marginals
    par(mfrow=c(2,2))
    for(i in 1:d) {
    for(j in 1:d) { #purposely not indented
        p = ks.test(sample$V[i,j,], ground$V[i,j,])$p.value #NB: order that ks.test takes its args is swapped from our order
        message(alg, " V[",i,",",j,"]", ": ", if(p > 0.05) { "same" } else { "different" })
    }
    }
}

#TODO: plot.NIW.marginals = plot.marginals
#TODO: write a function (to base plot.marginals on) which extracts the marginals in some sensible way;
#             either ruby-style where you pass an op block which receives the marginal as an argument
#             or more (procedurally?) where the marginals are reshaped into I guess a list with entries named by their marginal

#TODO: NIW.moment.first
# --> mean?
#TODO: NIW.moment.second
# --> ??

# run this all on NIW.naive
