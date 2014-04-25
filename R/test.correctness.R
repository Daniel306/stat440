
source("NIW.R")

#typedef: a multidimensional sample (MultiS) is
#  a matrix M with dim(M) = c(c(a,b, c, ...), n) 
#  n is the number of samples
#    (a, b, c, ...) is the dimensionality of the object

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
  if(is.function(ground)) {
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
  # if(!is.function(ground)) {
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



############
## Analytic densities
# we do not know analytic forms for many of the densities
# but we do know these two:

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

# TODO: NIW.X.marginal()

##########################################
## Moments


# utilties

cumulate <- function(M, margin, f) {
  # apply 'f' across 'margin'
  # this generalizes "cumsum", but so much slower.
  # f must take a vector and return a scalar
  # unlike apply(), this function protects the dimensions of M
  # this function is slow: it necessarily has at least quadratic runtime
  # use case: taking the cumulative variance in a (p,q,n) matrix (a matrix of n (p,q) design matrices, say)
  #  : cumulate(M, -3, var)
  # TODO: provide a feature so that you can flag to cumulate that f is already in cumulator form (eg so you can pass cumsum directly);
  #   this hope being that this will be faster
  
  original.dims = dim(M)
  M = apply(M, margin, function(v) { 
            sapply(1:length(v), function(i) { f(v[1:i]) })
            })
  dim(M) = original.dims #XXX if 
  return(M)
}



plot.convergence <- function(ground, samples, accumulator, title=NULL, ...) {
  #
  # note: plot.convergence does not actually guarantee you will see convergence; 
  #  you might not have enough samples or your choice of accumulator might not pick a duck out of a hat (an i.i.d. set of ducks, of course) and measure its beak length with some random error.
  #  using accumulator=mean or accumulator=var will definitely
  # ground and samples may be matrices
  # ... : args to plot
  # TODO: 

  # coerce ground to have the same shape as samples
  # ...
  # deal with the special case of ground being a single value (use case: an analytic reference value)
  # by giving it an extra dimension of length 1
  # TODO: factor this (how??) -- this step is like an inverse of drop(), except that drop() can't have a proper inverse
  if( length(dim(ground)) == length(dim(samples)) ) {
    # acceptable
  } else if(length(dim(ground)) == length(dim(samples)) - 1) { #TODO: also check the dimensions that match actually do match
    # fixable
    dim(ground) = c(dim(ground), 1)
  } else {
    # unacceptable!!
    stop("Inconsistent dimenesions between ground and sample; ground is c(", paste(dim(ground), sep=","),
                                                    "), while sample is c(", paste(dim(sample), sep=","), ")", sep="")
  }
  
  marginalized_dimension = length(dim(samples))

  # extract the name of the cumulation function for labelling
  accumulator.name = deparse(substitute(accumulator))
  
  if(missing(ylab)) {
    ylab = paste("sample", accumulator.name) #use  unless the user overrides
  }
  
  # Plot the means converging properlike
  kept_dimensions = (1:marginalized_dimension)[-marginalized_dimension] #precompute the list opposite of marginalized_dimension
                                             #doing it this way ensures this works even if ground is missing the last dimension
  samples = cumulate(samples, kept_dimensions, accumulator)
  ground = apply(ground, kept_dimensions, accumulator) # flattens
  dim(ground) = kept_dimensions
  stopifnot(dim(ground) == dim(samples)[1:2])
  
  marginals_do(samples, function(idx, samples) {
    plot(samples, xlab="samples taken (n)", ylab="sample mean", main=paste(title, idx2str(idx)))
    abline(h=do.call(`[`, c(ground, as.list(idx))), lty="dashed", col="blue")
    legend(paste(paste("True", accumulator.name), round(ground, 2)), lty="solid", col="blue")
  })
}


######################
## Computational (that is, kernel-density-estimated) Normal|Inverse-Wishart marginal moments 

plot.NIW.moment.first.computational <- function(ground, sample) { 
  # given a 'ground' sample (eg. from the naive implementation), compare
  # the marginals of the first (matrix-)moments of the rNIW samples visually
  plot.moment.first(ground$X, sample$X, "NIW X")
  plot.moment.first(ground$V, sample$V, "NIW V")
}


plot.NIW.moment.second.computational <- function(ground, sample) {
  # given a 'ground' sample (eg. from the naive implementation), compare
  # the marginals of the second (matrix-)moments of the rNIW samples visually
  plot.moment.second(ground$X, sample$X, "NIW X")
  plot.moment.second(ground$V, sample$V, "NIW V")
}

###########
#### Analytic Normal|Inverse-Wishart marginal moments 

 
NIW.mean <- function(Mu, Kappa, Psi, df, samples) {
  # compute the expected true mean of a 
  #  (X,V) ~ NIW(Mu, Kappa, Psi, df) distribution
  # 
  # returns: a list containing
  #   X -- the mean value for X
  #   V -- the mean value for V
  d = dim(Psi)[1]
  
  # the means of all Xs is just Mu, because X | V ~ N(Mu, V)
  # and the means of all Vs is Psi scaled by its degrees of freedom (props to Wikipedia for this one)
  list(X = Mu, V = Psi/(df-p-1))
}


NIW.var <- function(Mu, Kappa, Psi, df, samples) {
  # plot
  # compute the expected true variance
  #  (X,V) ~ NIW(Mu, Kappa, Psi, df) distribution
  #
  # does NOT return covariances; for one thing, some of those covariances have never been derived
  #  for another: no one cares
  # 
  # returns: a list containing
  #   X -- the variances of X
  #   V -- the variances of V
  
  d = dim(Psi)[1]

  # the Xs are t-distributed
  # which means their variances are ...
  # ..TODO
  X_var = NULL
  
  # the variances of the entries are
  # Var(V_ij) = (df-d+1)Psi_ij^2 + (v-p-1)*Psi_ii*Psi_jj) / (df-d)(df-d-1)^2(df-d-3)
  #
  # TODO: vectorize this formula
  # part of it vectorizes nicely, but the way to do the rest 
  # is not jumping out at me.
  # so a loop it is
   # cite: https://en.wikipedia.org/wiki/Inverse-Wishart_distribution
  V_var = (df- d + 1)*Psi^2
  for(i in 1:d) {
  for(j in 1:d) {
    V_var[i,j] = V_var[i,j] + (df - p - 1) * Psi[i,i]*Psi[j,j]
  }}  
  V_var = V_var / (df-d) / (df-d-1)^2 / (df-d-3)
  
  list(X = X_var, V = V_var)
}



##########################################
## Kolmogorov-Smirnov Tests

NIW.ks.test <- function(ground, sample) {
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
  
  # X marginals
  par(mfrow=c(2,2))
  for(i in 1:d) {
    p = ks.test(sample$X[i,], ground$X[i,])$p.value #NB: order that ks.test takes its args is swapped from our order
    message("\t X[",i,"]", ": ", if(p > 0.05) { "same" } else { paste("different", " (", p, ")", sep="")  })
  }
  
  # V marginals
  par(mfrow=c(2,2))
  for(i in 1:d) {
    for(j in i:d) { #purposely not indented
      p = ks.test(sample$V[i,j,], ground$V[i,j,])$p.value #NB: order that ks.test takes its args is swapped from our order
      message("\t V[",i,",",j,"]", ": ", if(p > 0.05) { "same" } else { paste("different", " (", p, ")", sep="") })
    }
  }
}
