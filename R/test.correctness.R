
source("NIW.R")

#typedef: a multidimensional sample (MultiS) is
#  a matrix M with dim(M) = c(c(a,b, c, ...), n) 
#  n is the number of samples
#    (a, b, c, ...) is the dimensionality of the object

# TODO: figure out a bloody better way of labelling the various output with which algorithm created it.


# utilties

cumulate <- function(s) {
  # wrap a statistic s(x) into a cumulative version
  #
  # this generalizes "cumsum", but so much slower.
  #
  # this function is slow: it necessarily has at least quadratic runtime,
  #  maybe more depending on how you write s.
  #
  # Do you understand that this function is slow, yet?
  #
  # args:
  #  s: must take a vector and return a scalar

  # use case: taking the cumulative variance in a (p,q,n) matrix (a matrix of n (p,q) design matrices, say)
  #  : apply(M, -3, cumulate(var))
  # if you can, use a more efficient specialized method
  
  function(v) { 
    sapply(1:length(v), function(i) { s(v[1:i]) })
  }
}

marginal.title <- function(title, idx) {
  # given one of my sketchy marginal indexing vectors
  # and a title prefix (which may be null)
  # construct a label that can be used to communicate what marginal we're displaying
  # SINCE this function is meant to be used in conjunction with marginals_do,
  # the last element of idx MUST be NA
  
  stopifnot(is.na(idx[length(idx)]))
  idx = idx[-length(idx)] # but we don't want to print that element
  
  paste(title, idx2str(idx), sep="") #sep="" means that a NULL title is a no-op
}


##################################################
## Density Plots

plot.density <- function(ground, sample, ...) { #XXX NAME
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

plot.densities <- function(ground, sample, title=NULL, layout=c(2,2), ...) { #XXX name
  # compare the marginals of the given sample to their expected p.d.f. curves
  #
  # this is basically just a loop over plot.density(),
  #  but it usefully accounts for labelling the plots uniformly
  #
  # args:
  #   ground: either, a MultiS, the same dimensionality as sample, in which case
  #         ground is assumed to be samples from the true distribution, and creates kernel density estimates
  #       -or-
  #       (( some kind of weird matrix-list containing p.d.f.s )) 
  #   sample: a MultiS or a vector
  #   layout: how to pack the plots (row by column)
  #   title: optional title to label thigns with
  #       
  
  # TODO: typechecks
  par(mfrow=layout)
  marginals_do(sample, function(idx, sample) {
    # case 1: ground is a matrix of samples
    plot.density(..getitem..(ground, idx), sample, main=marginal.title(title, idx), ...)
    # case 2:
    # TODO: ground is a matrix of pdfs
  })   
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
  list(X = Mu, V = Psi/(df-d-1))
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
  # Var(V_ij) = (df-d+1)Psi_ij^2 + (v-d-1)*Psi_ii*Psi_jj) / (df-d)(df-d-1)^2(df-d-3)
  #
  # TODO: vectorize this formula
  # part of it vectorizes nicely, but the way to do the rest 
  # is not jumping out at me.
  # so a loop it is
   # cite: https://en.wikipedia.org/wiki/Inverse-Wishart_distribution
  V_var = (df- d + 1)*Psi^2
  for(i in 1:d) {
  for(j in 1:d) {
    V_var[i,j] = V_var[i,j] + (df - d - 1) * Psi[i,i]*Psi[j,j]
  }}  
  V_var = V_var / (df-d) / (df-d-1)^2 / (df-d-3)
  
  list(X = X_var, V = V_var)
}



plot.convergence <- function(ground, samples, statistic, accumulator, title=NULL, ylab=NULL, ...) {

  # ground and samples may (should?) be matrices
  #  statistic is applied across the marginals of ground
  #  AS A SPECIAL USEFUL EXCEPTION, if ground only has one sample (ie it's a vector or a (d,1) matrix) 
  #   then and only then, ground is treated as the given expected value for statistic
  # AS ANOTHER SPECIAL EXCEPTION: if ground is NULL, no blue line is printed
  # note: plot.convergence does not actually guarantee you will see convergence; 
  #  you might not have enough samples or your choice of accumulator might not pick a duck out of a hat (an i.i.d. set of ducks, of course) and measure its beak length with some random error.
  #  using accumulator=mean or accumulator=var will definitely
  # IT IS UP TO YOU TO MAKE SURE statistic AND accumulator ARE COMPUTING THE SAME STATISTIC.
  
  #
  # ... : args to plot()

  #message("plot.convergence") #DEBUG
  if(is.null(ground)) {
    warning("ground is null; what are you doing with your life? the blue lines will not be plotted for comparison")
  }
  
  if(is.null(samples)) {
    stop("samples is null")
  }

  stopifnot(is.function(statistic))
  stopifnot(is.function(accumulator))

  # extract the name of the cumulation function for labelling
  statistic.name = deparse(substitute(statistic))
  
  if(is.null(ylab)) {
     ylab = paste("sample", statistic.name) #use  unless the user overrides
  }
  
  # coerce ground to have the same shape as samples
  # ...
  # deal with the special case of ground being a single value (use case: an analytic reference value)
  # by giving it an extra dimension of length 1
  # TODO: factor this (how??) -- this step is like an inverse of drop(), except that drop() can't have a proper inverse
  if(!is.null(ground)) {
     # TODO: INDENT THIS
  if( length(dim(ground)) == length(dim(samples)) ) {
    # acceptable
  } else if(is.null(dim(ground))) {
    # we have a vector; also fixable
    dim(ground) = c(length(ground), 1)
  } else if(length(dim(ground)) == length(dim(samples)) - 1) { #TODO: also check the dimensions that match actually do match
    # fixable
    dim(ground) = c(dim(ground), 1)
  } else {
    # unacceptable!!
    stop("Inconsistent dimenesions between ground and sample; ground is c(", paste(dim(ground), collapse=","),
                                                    "), while sample is c(", paste(dim(samples), collapse=","), ")", sep="")
  }
  }
  
  marginalized_dimension = length(dim(samples))
  
  # Plot the means converging properlike
  kept_dimensions = (1:marginalized_dimension)[-marginalized_dimension] #precompute the list opposite of marginalized_dimension
                                             #doing it this way ensures this works even if ground is missing the last dimension
  
  # Typecheck
  stopifnot(is.null(ground) || dim(ground)[kept_dimensions] == dim(samples)[kept_dimensions])

  #message("Looping") #DEBUG  
  marginals_do(samples, function(idx, samples) {
    samples = accumulator(samples)
    
    if(any(is.na(samples))) { return } #variances have their first element being NA, which causes problems; this check avoids all such problems with a hammer against a tree trunk
    
    # plot the converging samples against a reference horizontal line
    # and print the location of the horizontal line in the legend
    samples = reasonable_subset(cbind(1:length(samples), samples), 75) #plotting all the points slows down rendering and isn't enlightening;
                                    # Note:! just because we don't use all the points doesn't mean we can avoid the whole accumulator() call, becausethe last point is dependent on all the previous points
                                    # XXX TODO: if we're only going to print 75 points anyway, maybe it isn't unreasonable to use the quadratic algorithm *at each of those specific points*
                                    # and then we can make the API nicer and kill plot.*.convergence()
    plot(samples, xlab="samples taken (n)", ylab=ylab, main=marginal.title(title, idx))        
    
    if(!is.null(ground)) { # special case: ground can be NULL!
      ground = ..getitem..(ground, idx)
      if(length(ground) > 1) {   #special case: DON'T stat the ground truth if there's only one of it:
                                 #instead, assume it is the desired value of that statistic
        ground = statistic(ground)
      }
      stopifnot(length(ground) == 1) # at this point, ground MUST be scalar
      abline(h=ground, lty="solid", col="blue")
      legend("topleft", paste("True", statistic.name, "=", round(ground, 2)), lty="solid", col="blue")
    }
  })
}


plot.mean.convergence <- function(ground, samples, title=NULL, ...) {
  plot.convergence(ground, samples, mean, cummean, title=paste("Sample mean of", title), ...);
}

plot.var.convergence <- function(ground, samples, title=NULL, ...) {
  plot.convergence(ground, samples, var, cumvar, title=paste("Sample variance of", title), ...);
}



##########################################
## Kolmogorov-Smirnov Tests

ks.tests <- function(ground, sample, title=NULL, alpha=0.05) {
  # ground and sample must both be MultiSs of the same dimensionality
  # title is as in plot.densities and plot.convergence
  # alpha is the p-value cutoff for "different"
  marginals_do(sample, function(idx, sample) {
    ground = ..getitem..(ground, idx)   #only overwrites locally
    p = ks.test(sample, ground)$p.value #NB: ks.test() takes its args swapped from our convention
    message("\t", marginal.title(title, idx), ": ", if(p > alpha) { "same" } else { paste("different", " (", p, ")", sep="")  })
  })
}


NIW.ks.test <- function(ground, sample) {

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
  
  ks.tests(sample$X, sample$X, title="X")
  ks.tests(sample$V, sample$V, title="V")
}
