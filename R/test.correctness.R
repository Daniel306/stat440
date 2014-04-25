
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




############
## Analytic densities
# we do not know analytic forms for many of the densities
# but we do know these two:

IW.marginal.density <- function(i, j, Psi, df) {
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
    warning("non-diagonal marginals not supported (EVER)")
    NA
  }
}

dt.density <- function(df, mu, sigma.2) {
  # produce particular t-distribution density function
  # this is a simple partial evaluation wrapper around dt() from base
  #(I am worried the closure won't work right unless this is a whole separate function
  sd = sqrt(sigma.2)
  #message("dt.density closing over (df,mu,sd) = ", df, ", ", mu, ",", sd) #DEBUG
  function(x) {
    dt((x - mu) / sd, df)
  }
}

NIW.densities <- function(Mu, kappa, Psi, df, samples) {
  # [BROKEN]
  #  - $V is alright
  #  - $X is wrong; the variances come out screwy
  #  and it is weird that we are passing samples into an analytic density function -- shouldn't it just know from the params?
  #  (this is a little like saying to design a binary search tree you first need a running binary search tree)
  #  I do not know what the 3-param notation means; wikipedia only knows about up to a 2-param t dist'n: https://en.wikipedia.org/wiki/Noncentral_t-distribution
  #  in fact, isn't the point of the t that you *don't have to care* about the scale? that the variance cancels out?
  #  ----> but that's only true for .. hm.
  
  # n is the number of samples taken(?)
  # t(X) is the samples taken
  # it matters because the distribution of samples changes as more samples are taken
  # it only matters to the X side of the sampling, though, for the V side is i.i.d.

  X = samples$X;
  V = samples$V;
  X.bar = rowMeans(X)
  #message("X.bar");  print(dim(X.bar)); print(X.bar)

  #message("alloc ing the scatter matrix")
  # the Scatter matrix
  #print(dim(X))
  S = tcrossprod(X) #<-- the transpose is necessary because our convention (d,n) is opposite the usual
  #message("Scatter matrix")
  #print(dim(S))
  S.hat = tcrossprod(X.bar - Mu) # probably the wrong name for this quantity
  #print(S)
  #message("the othjer scatter matrix")
  #print(S.hat)
  
  d = dim(Psi)[1]
  n = dim(X)[1]
  
  X_pdf = as.list(rep(NA, d))
  V_pdf = as.list(rep(NA, d*d))
  dim(X_pdf) = c(d,1)
  dim(V_pdf) = c(d,d)

  # the marginals of the iWish part are iGamma on the diagonal
  # unknown elsewhere
  for(i in 1:d) {
    V_pdf[[i,i]] = IW.marginal.density(i,i, Psi, df)
  }

  for(i in 1:d) {
    df.t = df+n-d+1
    #y.bar = 0   # this is on the formula sheet by Kevin Murphy but not defined.. he probably meant x.bar.
                 # but even still, that is frustrating: it requires us to pass the data matrix in
    mu.t = (kappa*Mu[i] + n*X.bar[i])/(kappa+n)
    var.t = (Psi[i,i] + S[i,i] + kappa*n*S.hat[i,i]/(kappa+n))  /kappa/(df.t)
    X_pdf[[i, 1]] = dt.density(df.t, mu.t, var.t)
  }
  
  list(X = X_pdf, V = V_pdf)
}


plot.density <- function(ground, sample, breaks=40, ...) { #XXX NAME
  # plot histograms of the marginals in sample, overlaid with "ground truth" probability density.
  # 
  # ground: the "ground truth"
  #  this can either be a vector of samples, which will be kernel-density estimated, or a pdf function (which must be vectorized!!)
  #  the former is quicker to use: just generate samples in a naive and easy to validate way
  #  But the latter cannot be tainted by mistakes in your naive generator.
  # sample: the "sample"
  #  this is is the sample
  # breaks: 'breaks' parameter of hist() (see help(hist) for usage)
  # ...: extra arguments to plot()

  stopifnot(is.vector(sample))
  
  n = length(sample)
  
  # display the histogram
  hist(sample, probability=TRUE, breaks=breaks, xlab="", ...)

  if(is.null(ground) || is.na(ground)) {
    warning("plot.density got a null ground; not printing the pdf overlay")
    return(); # bail
  }
  
  # overlay the pdf
  if(is.function(ground)) {
    sample = reasonable_subset(sample, 77) #Reduce the number of points plotted
    sample = sort(sample)                  # this set must be sorted so that lines() dont zigzag everywhere
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
    #message("loop @ idx=") #DEBUG
    #print(idx)            #DEBUG
    title = marginal.title(title, idx) #Hack: do this first, because we possibly edit idx next
    
    if(length(dim(ground)) != length(idx)) { #HACK: assume this can only happen in the one-sample-special-casemagictime
 	                              # (this is a sketchy assumption)
      idx = rev(rev(idx)[-1]) #drop the last dimension because it's 
    }
    stopifnot(length(dim(ground)) == length(idx))
    ground = ..getitem..(ground, idx)
    if(class(ground) == "list") { #special case: you can pass a 'matrix' of density functions, by passing a list 
      stopifnot(length(ground) == 1)
      ground = ground[[1]]
      if(is.na(ground)) {
        ground = NULL
      } else {
        stopifnot(class(ground) == "function")
      }
    }
    plot.density(ground, sample, main=title, ...)
  })   
}


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



plot.moments <- function(ground, samples, statistic, accumulator, title=NULL, ylab=NULL, layout=c(2,2), ...) {

  # ground and samples may (should?) be matrices
  #  statistic is applied across the marginals of ground
  #  AS A SPECIAL USEFUL EXCEPTION, if ground only has one sample (ie it's a vector or a (d,1) matrix) 
  #   then and only then, ground is treated as the given expected value for statistic
  # AS ANOTHER SPECIAL EXCEPTION: if ground is NULL, no blue line is printed
  # note: plot.moments does not actually guarantee you will see moments or convergence to them
  #  you might not have enough samples or your choice of accumulator might not pick a duck out of a hat (an i.i.d. set of ducks, of course) and measure its beak length with some random error.
  #  using accumulator=mean or accumulator=var will definitely
  # IT IS UP TO YOU TO MAKE SURE statistic AND accumulator ARE COMPUTING THE SAME STATISTIC.
  #  and that that statistic is a moment
  
  #
  # ... : args to plot()

  par(mfrow=layout)

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
    plot(samples, xlab="samples taken (n)", ylab=ylab, main=marginal.title(title, idx), ...)
    
    if(is.null(ground)) { return() } # special case: ground can be NULL!
    
    ground = ..getitem..(ground, idx)
    if(length(ground) > 1) {   #special case: DON'T stat the ground truth if there's only one of it:
                                 #instead, assume it is the desired value of that statistic
      ground = statistic(ground)
    }
    stopifnot(length(ground) == 1) # at this point, ground MUST be scalar
    
    if(!is.finite(ground)) { return() } # bail if the number if bad at this point

    # finally, if we actually have a scalar number to use, use it!
    # TODO: un-hardcode the lty= and col=
    abline(h=ground, lty="solid", col="blue")
    legend("topleft", paste("True", statistic.name, "=", round(ground, 2)), lty="solid", col="blue")
  })
}


plot.means <- function(ground, samples, title=NULL, ...) {
  plot.moments(ground, samples, mean, cummean, title=paste("Sample mean of", title), ...);
}

plot.vars <- function(ground, samples, title=NULL, ...) {
  plot.moments(ground, samples, var, cumvar, title=paste("Sample variance of", title), ...);
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
