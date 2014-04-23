
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


#TODO: rewrite plot.NIW.marginals to it can just be 'plot.marginals'
#TODO: ditto for the ks tests
#  --> the trick used with moments will help here, where we write a generic function
#     which takes the fiddly details, main=, sub=, title=, etc, as arguments to be set by the caller
#    NIW.ks.test, in particular, definitely does not need to know what algorithm generated the data it is looking at
#    and that is a holdover from that we passed that argument in to plot.NIW.densities so that it could label the plots (because there's no other way in R but to label the plots as you create them)
#TODO: write a function (to base plot.marginals on) which extracts the marginals of a matrix in some sensible way;
#             either ruby-style where you pass an op block which receives the marginal as an argument
#             or more (procedurally?) where the marginals are reshaped into I guess a list with entries named by their marginal
#  the core of this function now exists as R() inside of plot.converging.moment.multi; factoring it out will be helpful: it will drastically shorten *.marginals() above
# TODO: figure out a bloody better way of labelling the various output with which algorithm created it.

##################################################
## Density Plots

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
  
  # display the histogram
  hist(sample, probability=TRUE, breaks=40, ...)
  
  # overlay the pdf
  if(is.call(ground)) {
    lines(sample, ground(sample), col="blue", lty="twodash")
  } else {
    # if caller has not given us a pdf directly, then
    # assume ground is a vector for us to kernel-density estimate
    stopifnot(is.vector(ground) && is.numeric(ground)) #typecheck
    
    lines(density(ground), col="blue", lty="twodash")
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


plot.NIW.densities <- function(ground, sample, alg=NULL) { #XXX name
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
  # and we also know that V is symmetric, so we can do only (1:d)x(i:d)
  # --> this tidbit will be difficult to factor out
  
  # X marginals
  par(mfrow=c(2,2))
  for(i in 1:d) {
    plot.compare(ground$X[i,], sample$X[i,], main=paste("X[",i,"]"), sub=alg)
  }
  
  # V marginals
  par(mfrow=c(2,2))
  for(i in 1:d) {
    for(j in i:d) { #purposely not indented
      plot.compare(ground$V[i,j,], sample$V[i,j,], main=paste("V[",i,",",j,"]"), sub=alg)
    }
  }
}



IW.marginal <- function(i, j, Psi, df) {
  # Inverse Wishart marginals
  #
  # returns:
  #   a vectorized function in 'x' giving the pdf giving the inverse wishart
  stopifnot(is.symmetric(Psi))
  d = nrow(Psi)
  
  if(i == j) {
    # from wikipedia
    # Any subset of rows/cols our things are Inverse Wishart:
    # SubW ~ IW(SubPsi, df - (number of columns outside of))
    # which isn't all that useful for us interesting.
    # but it does mean that a particular, scalar, diagonal entry
    # is
  
    # close over the scalar invgamma parameters instead of over a whole matrix
    # ((though maybe R closes over Psi anyway))
    
    alpha = (df - d - 1)/2
    beta  = (Psi[i , i])/2
    function(x) { dinvgamma(x, shape=alpha, scale=beta) }
  } else {
    stop("non-diagonal marginals not supported yet")
  }
}


##########################################
## Moments

# Wikipedia gives some useful analytic formulas:
# https://en.wikipedia.org/wiki/Inverse-Wishart_distribution#Moments
# and it's relatively obvious that the mean of X is just Mu.

plot.converging.moment <- function(ground, samples, ...) {
  # 
  # every moment is an expectation, which can be approximated by a sample mean
  # to use this function, first compute
  # plots a horizontal line at the expected value as given in 'ground
  #
  # args:
  #  ground: the 'ground truth': a number or a vector
  #  samples: a vector of sample points
  #  ...: extra arguments to plot()

  # TODO: typechecks
  
  if(length(ground) > 1) { ground = mean(ground) }

  # plotting too many points causes lag
  # so use intseq() to reduce
  idx = intseq(1, length(samples), length.out=1000)
  means = cummean(samples)[idx]
  plot(idx, means, xlab="samples taken", ylab="moment", ...) #XXX is ylab right?
  abline(h=ground, lty="solid", col="blue")
}


plot.converging.moment.multi <- function(ground, samples, title=NULL) { # XXX name
  # befpore using this function, map samples (and ground) with
  # for example, to take second moments, first run crossprod() across samples
  #
  # TODO: docstring
  # TODO: merge plot.converging.moment into this, so that it handles multivariates smoothly in line with univariates?
  #    con: the plotting arguments, main=, sub=, etc, will be finicky to do if it's merged
  #    pro: it might (might) be faster to do the cummeans all at once

  # TODO: typechecks
  
  # our convention is: the LAST dimension is 'n'  
  n_dimension = length(dim(samples))
  interesting_dimensions = 1:(n_dimension-1) #XXX name
  
  #ground = apply(ground, -n_dimension, mean) # ground is FLATTENED to a single number (per marginal, so there's actually still a bunch)
  #dim(ground) = dim(samples)[interesting_dimensions] #for some reason, using mean() means the result ground loses all its dimensions and becomes a vector. go figure.
  #samples = apply(samples, -n_dimension, cummean) # PAY ATTENTION: this line changes samples from (d,d,n) to (n,d,d)
    
  # deal with the special case of ground being a single value
  # by giving it an extra dimension of length 1
  # (and the beginnings of typechecks while we're at it) 
  # ...TODO: should this part be killed? It's replicating work that happens in functions that call it.
  if( length(dim(ground)) == length(dim(samples)) ) {
    # acceptable
  } else if(length(dim(ground)) == length(dim(samples)) - 1) { #TODO: also check the dimensions that match actually do match
    # fixable
    dim(ground) = c(dim(ground), 1)
  } else {
    # unacceptable!!
    stop("Inconsistent dimenesions between ground and sample")
  }

  # loop over all interesting_dimensions and plot each marginal
  # recursive something someting  
  R<-function(idx, dims) { #TODO: pull out into util.R
    if(length(dims)==0) {
      # base case
      # we've bottomed out and have an actual index in hand
      # so use it
      #message("bottomed out at idx=")
      #print(idx) #DEBUG
      
      # FIXME: I don't understand how to handle variable arity functions in R
      #  (I would really like to say "samples[*idx]" like I can in python)
      # so I've hardcoded the cases we're actually using
      # this needs to be repaired!
      if(length(idx) == 0) {
        plot.converging.moment(ground[], samples[],  main=paste(title, sep=""));
      } else if(length(idx) == 1) {
        i = idx[1];
        plot.converging.moment(ground[i,], samples[i,],  main=paste(title,"[",i,"]", sep=""));
      } else if(length(idx) == 2) {
        i = idx[1]; j = idx[2];
        plot.converging.moment(ground[i,j,], samples[i,j,],  main=paste(title, "[",i,",",j,"]", sep=""));
      }
    } else {
      # split the top dimensions from the body
      # note how we extract the *number* of elements in the dimension d from the index of the dimension itself dims[1]
      d = dim(samples)[dims[1]]; dims = dims[-1];
      for(i in 1:d) {
        R(append(idx, i), dims)
      }
    }
  }
  R(c(), interesting_dimensions)
}
# sketchy test:
#test.plot.converging.moment.multi <- function() {
#  source("test.constants.R")
#  S = rWishart(14441, kDF, kPsi);
#   # the mean of a wishart distribution is V*df..
#  plot.converging.moment.multi(kPsi*kDF, S, paste("W(Psi,", kDF,")"))
#}
#test.plot.converging.moment.multi()


moment.first <- identity


moment.second <- function(M) {
  # given a (d,n) or (d,p,n) matrix
  #  (n being the number of samples)
  # compute the second moment of each sample
  # which, in matrixland, is M%*%t(M), a (d,d) matrix (nb: a (d,n) matrix means each sample is (d,) which tcrossprod treats as a (d,1) matrix)
  # returns: a (d,d,n) matrix containing the second moments
  # TODO: generalize to matrices of more than 2 dimensions

  # Typechecks
  n_dim = length(dim(M))
  d = dim(M)[1]
  n = dim(M)[n_dim]
  
  R = apply(M, n_dim, tcrossprod)
  # apply() flattens its results just to be a pain;
  # we need to unflatten them.
  # M %*% t(M) is (d x p)*(p * d) so the result is (d x d)
  dim(R) = c(d, d, n)
  
  return(R)
}



plot.moment.first <- function(ground, samples, title=NULL) {
  if(length(dim(ground)) == length(dim(samples)) - 1) {
    dim(ground) = c(dim(ground), 1)
  } # else: maybe unhandled crashytimes
  
  plot.converging.moment.multi(moment.first(ground), moment.first(samples), title=title)
}


plot.moment.second <- function(ground, samples, title=NULL) {
  #follow along with how finicky this is:
  # the mainline case is that ground and sample are both (d,p,n) matrices
  # but ground MIGHT be a single sample, without the 'n' dimension
  # we special- this case by giving it a third dimension of length 1
  
  if(length(dim(ground)) == length(dim(samples)) - 1) {
    dim(ground) = c(dim(ground), 1)
  } # else: maybe unhandled crashytimes
  
  plot.converging.moment.multi(moment.second(ground), moment.second(samples), title=paste(title, "^2", sep=""))
}


plot.NIW.moment.first <- function(ground, sample) { 
  #TODO: should this take the distribution paramters so it can compute 'ground' directly?
  plot.moment.first(ground$X, sample$X, "NIW X")
  plot.moment.first(ground$V, sample$V, "NIW V")
}


plot.NIW.moment.second <- function(ground, sample) {
  #TODO: should this take the distribution paramters so it can compute 'ground' directly?
  plot.moment.second(ground$X, sample$X, "NIW X")
  plot.moment.second(ground$V, sample$V, "NIW V")
}

# TODO: the formulas on wikipedia include the covariances of each element of V; we aren't computing covariances here, but maybe we should.
# TODO: some of the plots disagree with their ground truth. curious. the ks tests don't complain, though.

##########################################
## Kolmogorov-Smirnov Tests

NIW.ks.test <- function(ground, sample, alg=NULL) {
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
  # and we also know that V is symmetric, so we can do only (1:d)x(i:d)
  # --> this tidbit will be difficult to factor out
  
  message("KS Tests")
  # X marginals
  par(mfrow=c(2,2))
  for(i in 1:d) {
    p = ks.test(sample$X[i,], ground$X[i,])$p.value #NB: order that ks.test takes its args is swapped from our order
    message(alg," X[",i,"]", ": ", if(p > 0.05) { "same" } else { paste("different", " (", p, ")", sep="")  })
  }
  
  # V marginals
  par(mfrow=c(2,2))
  for(i in 1:d) {
    for(j in i:d) { #purposely not indented
      p = ks.test(sample$V[i,j,], ground$V[i,j,])$p.value #NB: order that ks.test takes its args is swapped from our order
      message(alg, " V[",i,",",j,"]", ": ", if(p > 0.05) { "same" } else { paste("different", " (", p, ")", sep="") })
    }
  }
}
