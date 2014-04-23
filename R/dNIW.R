# dNIW.R
#
# This file implements dNIW(), which gives the density function (p.d.f.)
#  of the Normal|Inverse-Wishart distribution.
#
# It also implements dMNIW(), which is the density function of the Matrix-Normal|Inverse-Wishart.
#

source("util.R")
source("NIW.util.R")

dNIW.typecheck <- function(X, V, Mu, Kappa, Psi, df) {
  NIW.typecheck(Mu, Kappa, Psi, df)
  
  # additionally, X and V must match the dimens
  # the above enforces that Mu is a vector
  
  stopifnot(is.vector(X) && length(X) == length(Mu))
  stopifnot(is.symmetric(V) && dim(V) == dim(Psi))
}


dNIW.typechecked <- function(dNIW) {
    # wrap an rNIW implementation "rNIW" in the common safety code that
    #  all rNIW implementations should share.
    # this enforces all the pre (and post) conditions necessary.
    # in particular, inside of rNIW it is safe to say "d = length(Mu)" (XXX would it be nice to pass d INTO rNIW...??)  
  
    function(X, V, Mu, Kappa, Psi, df, log=FALSE) {
       #PRECONDITIONS
        dNIW.typecheck(Mu, Kappa, Psi, df)
        
        # IMPLEMENTATION
        ans = dNIW(X, V, Mu, Kappa, Psi, df, log=log);
        
        # POSTCONDITIONS
        
        stopifnot(class(ans) == "numeric")
        stopifnot(length(ans) == 1)
        # enforce that probability result is in [0,1]
        # XXX is this correct? it *is* possible for probabilities to be larger than 1 in weird corner cases of continuousland, isn't it?
        if(log) {
          stopifnot(ans < 0)
        } else {
          stopifnot(0<=ans && ans<=1)
        }
        
        return(ans)
    }
}



dIW <- function(V, Psi, df, log=FALSE) {
  # TODO: vectorize this

  # Typechecks
    stopifnot(is.symmetric(Psi))  
      d = dim(Psi)[1]
        stopifnot(dim(V)[1:2] == c(d,d))

  # coerce V to a vectorized case, even when it isn't
    # this lets us handle everything uniformly
      if(length(dim(V)) == 2) {
          dim(V) = c(d,d,1)  
            }


  c = v*log(det(Psi))/2 - (df*d)*log(2)/2 - multigamma(d)(df/2) #normalizing constant

  det_V = apply(V, 3, det)
    #Psi %*% solve(V) would be a clever backsolve IF we had the square root of V handed to us
      trace_Psi_invV = apply(Psi_invV, 3, function(M) { sum(diag( Psi %*%  solve(M)    )) } )

  # Wikipedia is wrong as of <https://en.wikipedia.org/w/index.php?title=Inverse-Wishart_distribution&oldid=604297660>
    # the pdf on https://en.wikipedia.org/wiki/Inverse-Wishart_distribution should be identical to https://en.wikipedia.org/wiki/Wishart_distribution
      #  but with the X and V terms inverted
        # and indeed those are
          # but the Inverse Wishart page additionally switches the 2/3 signs in the power
            # which is completely the sort of mistake that one might make with LaTeX
              # ...but two software packages (LaplacesDemon and MCMCpack) agree with the wikipedia formula. so... hm.
                p = -  (   (df+d+1)*log(det_V)  + trace_Psi_invV    )/2    + c

  if(!log) {
      p = exp(p)
        }

  return(p)
  }

dNIW <- function(X, V, Mu, Kappa, Psi, df, log=FALSE) {
  # dNIW(X, V, Mu, Kappa, Psi, df, log=FALSE):
  #  produces probability densities P(NIW(Mu, Kappa, Psi, df) = (X, V))
  #
  # args:
  #     X, V: the point(s) at which to select:
  #      either:
  #       X is a d-length vector or a (d,1) matrix
  #       V is a (d,d) matrix
  #      or:
  #       X is a (d,n) matrix
  #       V is a (d,d,n) matrix
  #     Mu, Kappa, Psi, df: the parameters of the particular Normal/Inverse-Wishart distribution in question
  #     log: whether to produce log-probabilities or not.
  #
  # preconditions:
  #   XXX FILL ME IN
  #
  # returns:
  #   a vector of probabilities p, or log(p) if log==TRUE of length n (possibly n=1)

  #
    # Typechecks
      #TODO: test that X and V can be vectorized..

  # coerce the arguments to vectors

  # always compute in logspace, because it's less prone to roundoff
    # map to normal space at the end if requested

  # XXX both these ops could be much sped up if we
    # were given the square root of V instead
      pN = dmvnorm(X, Mu, V/Kappa, log=T)
        pIW = dIW(V, Psi, df, log=T)

  p = pN + pIW; # P(X,V) = P(X|V)P(V) = exactly the piece distributions used to generate the NIW in the first place.

  if(!log) {
      p = exp(p)
        }

  return(p)
  }
  