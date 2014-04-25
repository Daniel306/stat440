
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
# TODO: its pretty obvious these functions are identical and can be factored
# but how is eluding me right now
# the trouble is that these intimately depend on the knowledge that cummean/cumvar return a vector the same size as the input
# and I don't see how to safely factorize out over that
# a day with the apply manpage might change my mind.

cummean.marginals = function(M) {
  # precondition: M's dim is (a,b,c,d,...,n); marginals are taken of that last dimension, n
  # postcondition: a matrix the same size as M, with each of those marginal vectors replaced by its cumulative mean
  original.dims = dim(M)
  n_dim = length(original.dims)
  M = apply(M, -n_dim, cummean)
  dim(M) = original.dims #unflatten what apply() flattened
  return(M)
}

cumvar.marginals = function(M) {
  # precondition: 
  # precondition: M's dim is (a,b,c,d,...,n); marginals are taken of that last dimension, n
  # postcondition: a matrix the same size as M, with each of those marginal vectors replaced by its cumulative variance
  original.dims = dim(M)
  n_dim = length(original.dims)
  M = apply(M, -n_dim, cumvar)
  dim(M) = original.dims #unflatten what apply() flattened
  return(M)
}

marginals_do <- function(M, f) { #TODO: pull out into util
  
  # precondition: M's dim is (a,b,c,d,...,n); marginals are taken of that last dimension, n    
  # precondition: f is f(idx, v) where idx will be a vector containing the length(dim)-1 indexes and v is the marginal vector
  # postcondition: nothing; you need to do every with side-effects in f
  # loop over all interesting_dimensions and plot each marginal
  # recursive something someting  
  R<-function(idx, dims) { #TODO: pull out into util.R
    if(length(dims)==0) {
      # base case
      # we've bottomed out and have an actual index in hand
      # so use it
      #message("bottomed out at idx=")
      #print(idx) #DEBUG
      
      v = do.call(`[`, c(M, as.list(idx))) # this line courtesy of mrflick in #R
      v = as.vector(v) # just in case (i think if idx == c(), which should correspond to M[], we will *not* get a vector out: R will preserve the dimensions (but only in that case (which shouldn't be possible, if the recursion is correct))
      f(idx,  v)
    } else {
      # recursive case
      # split the top dimensions from the body
      # note how we extract the *number* of elements in the dimension d from the index of the dimension itself dims[1]
      d = dim(samples)[dims[1]]; dims = dims[-1];
      for(i in 1:d) {
        R(append(idx, i), dims)
      }
    }
  }
  R(c(), interesting_dimensions) #kick off the recursive loop  
  return(NULL)
}


reasonable_subset(D, length.out=1000) {
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


# TODO: the following two functions should be factored
# but, as with the above, how to do so eludes me
#  because again, there's special-case behaviour in each; plot.convering.variance would, for example, behave erratically if it ended up calling var(scalar) which is NA

plot.converging.mean <- function(ground, samples, ...) {
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

  n = length(samples)
  idx = 1:n

  plot(reasonable_subset(cbind(idx, cummean(samples)), means, xlab="samples taken", ylab="mean", ...) # XXX is ylab what we want?
  #  note: ^this relies on plots() special-case behaviour of plot on a 2-column matrix/dataframe being interpreted the same as plot(x,y)
  abline(h=ground, lty="solid", col="blue")
  legend(paste("True mean", round(ground, 2)), lty="solid", col="blue")
}


plot.converging.variance <- function(ground, samples, ...) {
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
  
  if(length(ground) > 1) { ground = var(ground) }

  n = length(samples)
  idx = 1:n
  
  plot(reasonable_subset(cbind(idx, cumvar(samples)), means, xlab="samples taken", ylab="mean", ...) # XXX is ylab what we want?
  #  note: ^this relies on plots() special-case behaviour of plot on a 2-column matrix/dataframe being interpreted the same as plot(x,y)
  abline(h=ground, lty="solid", col="blue")
  legend(paste("True variance", round(ground, 2)), lty="solid", col="blue")
}


idx2title <- function(idx) {
  title = paste(idx, collapse=",") #could also do this with a complicated do.call, but paste() covers this case helpfully for us with 'collapse'
  if(length(idx) > 0) { #special case: no indices means all the indecies #XXX untested
    title = paste("[", title, "]") #but otherwise wrap the indecies in square brackets
  }
  return(title)
}
        

plot.converging.moment.multi <- function(ground, samples, title=NULL, sub=NULL) { # XXX name
  # befpore using this function, map samples (and ground) with
  # for example, to take second moments, first run crossprod() across samples
  #
  # TODO: docstring
  # TODO: merge plot.converging.moment into this, so that it handles multivariates smoothly in line with univariates?
  #    con: the plotting arguments, main=, sub=, etc, will be finicky to do if it's merged
  #    pro: it might (might) be faster to do the cummeans all at once
  # TODO: factor the bit that loops down dimensions to a subroutine
  
  # TODO: typechecks
  
  # our convention is: the LAST dimension is 'n'  
  n_dimension = length(dim(samples))
  interesting_dimensions = 1:(n_dimension-1) #XXX name
  
  #ground = apply(ground, -n_dimension, mean) # ground is FLATTENED to a single number (per marginal, so there's actually still a bunch)
  #dim(ground) = dim(samples)[interesting_dimensions] #for some reason, using mean() means the result ground loses all its dimensions and becomes a vector. go figure.
  #samples = apply(samples, -n_dimension, cummean) # PAY ATTENTION: this line changes samples from (d,d,n) to (n,d,d)
    
  # deal with the special case of ground being a single value (use case: an analytic reference value)
  # by giving it an extra dimension of length 1
  
  if( length(dim(ground)) == length(dim(samples)) ) {
    # acceptable
  } else if(length(dim(ground)) == length(dim(samples)) - 1) { #TODO: also check the dimensions that match actually do match
    # fixable
    dim(ground) = c(dim(ground), 1)
  } else {
    # unacceptable!!
    stop("Inconsistent dimenesions between ground and sample")
  }

  marginals_do(samples, function(idx, v) {    

    # awkward: marginals_do is looping down samples, but we really are *simultaneously* looping over ground
    #  maybe some sort of zip()-like construct would help
    # for now, copying the magic line from marginals_do and just reaching up one scope works:
    ground_marginal = do.call(`[`, c(ground, as.list(idx)))
    
    plot.converging.moment(ground_marginal, v, main=paste(title, idx2title(idx)), sub=sub);
     # ^ TODO: factor this
    })
      
}
# sketchy test:
#test.plot.converging.moment.multi <- function() {
#  source("test.constants.R")
#  S = rWishart(14441, kDF, kPsi);
#   # the mean of a wishart distribution is V*df..
#  plot.converging.moment.multi(kPsi*kDF, S, paste("W(Psi,", kDF,")"))
#}
#test.plot.converging.moment.multi()

#######
## Computational matrix-moments

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


# Note well: the marginals of the first moment are the same as the first moments of the marginals
# BUT, this is not true in general

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

 
plot.NIW.moment.mean.analytic <- function(Mu, Kappa, Psi, df, samples) {
  # plot (???) and then, suddenly, rabbits.
  # this is not called ".first" because its cousin below plots variances, not second moments
  
  d = dim(Psi)[1]
  
  # the means of all Xs is just Mu, because X | V ~ N(Mu, V)
  plot.converging.moment.multi(Mu, samples$X, title="NIW X", sub="(expected value from analytic formulas)")
  
  # the means of all Vs is Psi scaled by its degrees of freedom
  plot.converging.moment.multi(Psi/(df-p-1), samples$V, title="NIW V", sub="(expected value from analytic formulas)")
}


plot.NIW.moment.variance.analytic <- function(Mu, Kappa, Psi, df, samples) {
  # plot the second moments
  # this is not called ".first" because its cousin below plots variances, not second moments
  
  d = dim(Psi)[1]
  
  # the Xs are t-distributed
  # which means their variances are ...
  # ..TODO
  
  # the variances of the entries are
  # Var(V_ij) = (df-d+1)Psi_ij^2 + (v-p-1)*Psi_ii*Psi_jj) / (df-d)(df-d-1)^2(df-d-3)
  #
  # TODO: vectorize this formula
  # part of it vectorizes nicely, but the way to do the rest 
  # is not jumping out at me.
  # so a loop it is
  V_var = (df- d + 1)*Psi^2
  for(i in 1:d) {
  for(j in 1:d) {
    V_var[i,j] = V_var[i,j] + (df - p - 1) * Psi[i,i]*Psi[j,j]
  }}  
  V_var = V_var / (df-d) / (df-d-1)^2 / (df-d-3)
  
  # cite: https://en.wikipedia.org/wiki/Inverse-Wishart_distribution
  plot.converging.variance.multi(v_var, samples$V, title="NIW variance of V", sub="(expected value from analytic formulas)")
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
