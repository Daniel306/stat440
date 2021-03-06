# rNIW.R
#
# This file implements several versions of the function rNIW(), which generates
# random samples from the Normal|Inverse-Wishart distribution
# <https://en.wikipedia.org/wiki/Normal-inverse-Wishart_distribution>.
#
# There is also a generalization, rMNIW(), which samples the
# Matrix-Normal|Inverse-Wishart distribution, which is the same as the NIW,
# except that the X sample points are themselves rectangular matrices as
# given by the Matrix-Normal distribution <https://en.wikipedia.org/wiki/Matrix_normal_distribution>.
#
# This is prototype code. Eventually, one version will be selected, packaged, and uploaded to CRAN.
#
# There is a standard API. In lieu of repeating ourselves, take the following as contained in the docstring for each implementation.
# 
# rNIW(n, Mu, Kappa, Psi, df):
#  generate samples (X, V) ~ NIW(Mu, kappa, Psi, df)
#
# see `FIXME` for the full definition and purpose of the NIW distribution. 
# args:
#     n: the desired number of sample "points" (each "point" is a pair of a d-length vector and a dxd matrix)
#     Mu:  a vector giving the mean of X samples
#     Kappa: a scalar which adjusts the variances of X
#     Psi: a matrix, the mean of V samples
#     df: the "degrees of freedom" of the underlying Wishart distribution
#
# preconditions:
#       the dimensions match: d := length(Mu) == dim(Psi)[1] == dim(Psi)[2];
#       kappa is positive.
#       Psi is is positive-definite (=> square, symmetric, all positive eigenvalues, and chol(Psi) exists).
#       df need not be integer, but it must be larger than d-1.
#
# returns:
#   a list with elements
#     $X, the c(d, n) matrix of sample X values.
#     $V, the c(d, d, n) matrix of sample V values.
#  $X and $V have the same number of samples, and ans$X[,i] corresponds is a d-dimensional normal generated with variance ans$V[,,i]
#  each V[,,i] ~ inv(W(inv(Psi), df))
#  each X[,i] ~ N(Mu, V[,,i]/Kappa)
# NB: having n as the last dimension means that it is the 'widest' dimension:
#     each set of d or d^2 sample elements run consecutively.
#     This format is as returned by the built-in rWishart() function.


source("util.R")
source("NIW.util.R")
#require("slam") #TODO: investigate this package

require("Rcpp")
Rcpp::sourceCpp("rNIW.cpp")

##############################################################################
#                          Helper functions
##############################################################################

rNIW.alloc <- function(n, d) {
  # Construct an empty datastructure to store rNIW results.
  #
  # args: 
  #   n: the number of samples to make space for
  #   d: the dimensionality of the
  #
  # returns:
  #   a list matching the rNIW return value spec
  # 
  # Beware: *pre*allocating will give a performance hit!
  #        R is call-by-value, so passing a matrix in an argument does a memcpy()
  # Thus this is only useful if used in the same function where the elements are filled.
  
  X = matrix(NA, d*n);
  dim(X) = c(d, n);

  V = matrix(NA, d^2*n);
  dim(V) = c(d,d,n);
  
  list(X=X, V=V)
}

rNIW.typecheck <- function(n, Mu, Kappa, Psi, df) {
  # rNIW preconditions
  
  # n must be a natural number
  stopifnot(length(n) == 1 && n==round(n) && n > 0) # is.integer() is wrong for this; see its help

  # the distribution parameters must be sane
  NIW.typecheck(Mu, Kappa, Psi, df)
}


rNIW.typechecked <- function(rNIW) {
    # wrap an rNIW implementation "rNIW" in the common safety code that
    #  all rNIW implementations should share.
    # this enforces all the pre (and post) conditions necessary.
    # in particular, inside of rNIW it is safe to say "d = length(Mu)" (XXX would it be nice to pass d INTO rNIW...??)  

    function(n, Mu, Kappa, Psi, df) {
        rNIW.typecheck(n, Mu, Kappa, Psi, df)
        d = nrow(Psi)
        
        # IMPLEMENTATION
        ans = rNIW(n, Mu, Kappa, Psi, df);

        # POSTCONDITIONS

        # make sure the results exist
        stopifnot(class(ans) == "list")
        stopifnot(all(sort(names(ans)) == c("V", "X")))

        # check dimensionality of the results
        stopifnot(dim(ans$X) == c(d,n))
        stopifnot(dim(ans$V) == c(d,d,n))

        # XXX similarly, we don't check that Vs are positive definite (or even symmetric), though they should be

        return(ans)
    }
}


rMNIW.typecheck <- function(n, Mu, Kappa, Psi, df) {
  # rMNIW preconditions
  
  # n must be a natural number
  stopifnot(length(n) == 1 && n==round(n) && n > 0) # is.integer() is wrong for this; see its help
  
  # the distribution parameters must be sane
  MNIW.typecheck(Mu, Kappa, Psi, df)
}

rMNIW.typechecked <- function(rMNIW) {
  # wrap an rMNIW implementation "rNIW" in the common safety code that
  #  all rNIW implementations should share.
  # this enforces all the pre (and post) conditions necessary.
  # in particular, inside of rNIW it is safe to say "d = length(Mu)" (XXX would it be nice to pass d INTO rNIW...??)  
  
  function(n, Mu, Kappa, Psi, df) {
    rMNIW.typecheck(n, Mu, Kappa, Psi, df)
    d = nrow(Psi)
    q = nrow(Kappa)
    # IMPLEMENTATION
    ans = rMNIW(n, Mu, Kappa, Psi, df);
    
    # POSTCONDITIONS
    
    # make sure the results exist
    stopifnot(class(ans) == "list")
    stopifnot(all(sort(names(ans)) == c("V", "X")))
    
    # check dimensionality of the results
    stopifnot(dim(ans$X) == c(q,d,n))
    stopifnot(dim(ans$V) == c(d,d,n))
    
    # XXX similarly, we don't check that Vs are positive definite (or even symmetric), though they should be
    
    return(ans)
  }
}

##############################################################################
#                           Implementations
##############################################################################


##################################
### Naive Sampling

require(MASS);
rmvnorm <- MASS::mvrnorm;  #repair naming convention;
mvrnorm <- NULL;           #I don't know why I need the namespace in the first step,
                           #but without it this breaks under source(). :shrug:

rInvWishart <- function(n, df, Sigma) {
 # naive inverse wishart sampler 
 stopifnot(nrow(Sigma) == ncol(Sigma))
 stopifnot(df >= nrow(Sigma))
 W = rWishart(n, df, solve(Sigma));
 for(i in 1:n) {
   W[,,i] = solve(W[,,i]) 
 }
 W
}

rNIW.extremelynaive <- rNIW.typechecked(function(n, Mu, kappa, Psi, df) {
  # generate the n samples by n times doing: i) generate V, ii) generate X given V.
  # this "extremely naive" version is a 1-for-1 translation of the algorithm on wikipedia
  #  the "extreme" part is that it generates every sample one at a time in an R loop
  
  d = length(Mu)
  ans = rNIW.alloc(n,d)
  
  for(i in 1:n) {
     ans$V[,,i] = rInvWishart(1, df, Psi)[,,1];
     ans$X[,i] = rmvnorm(1, Mu, ans$V[,,i]/kappa);
  }
  
  ans
})

rNIW.naive <- rNIW.typechecked(function(n, Mu, kappa, Psi, df) {
  # generate the n samples
  # this calls rWishart once to generate the n wishart samples, then inverts them all (to make them inverse wishart samples)
  # The only effective difference between this and extremelynaive should be that this only inverts Psi once
  
  # because we use rWishart in one step in this implementation
  # it is the only one in which we do not use rNIW.alloc()
  d = length(Mu)
  ans = list()
  ans$X = matrix(NA, d, n)
  ans$V = rInvWishart(n, df, Psi);   
  
  # unfortunately, generating X must still be in a loop because each sample has a different distribution
  # and there's no way with rmvnorm to ask "give me n samples and here's the parameters for each"
  # ..is there?
  #
  for(i in 1:n) {
     ans$X[,i] = rmvnorm(1, Mu, ans$V[,,i]/kappa);
  }
  # the other option is to generate the multivariate normals ourselves (based on standard normals: rnorm(d*n, 0, 1))
  # and put everything together in a single matrix multiply or something
  # but to do that requires making a blockwise-diagonal matrix, which means lots and lots of zeros being multiplied for no reason
  # it is probably more computation to do that than to just use a loop
  # and the above has the advantage of calling a reliable standard library function to do the difficult distribution calculation for us
  
  ans
})


# naive algorithm, but done in C, for a more fair comparison
rNIW.extremelynaive.RcppEigen <- rNIW.typechecked(rNIW_extremelynaive_eigen)


##################################
###   Sampling by exploiting
### the Bartlett Decomposition.

# The key idea is that the canonical way to generate Wishart matrices
#  first gets their cholesky decomposition, i.e. a square root of the matrix.
#  this is called the "Bartlett Factor" for it is Bartlett's Theorem
# that says this method works.
# And to generate MultiNormals with a particular covariance matrix,
#  you scale a standard normal Z by a square root of that matrix.
#  (rmvnorm is definitely doing just that).
# So if you need to generate the joint Normal-(Inverse)-Wishart distribution, 
#  there should be some speed gains by avoiding the useless square-squareroot.
#
# Additionally, once we solidify a prototype, we're going to port it to C,
#  which will further speed it up.

# confused, scrap notes on this process:
# is it better to actually construct U than try to apply it after the fact?
#   --> if all we wanted was X, then no, but we also need V, which means we MUST construct U somewhere
# (notice: U is the factor that is necessary...)
# TODO: only use U instead of using A in one place and U in another
   # what order is it?
   # t(U)%*%U = W; but we want W^-1 = V, so we need V = U^-1 
     # --> construct U?
     # --> U is
     # gamma
     
     #
     # and here's where this algorithm starts to look silly: we EITHER invert gamma and then invert it back, or we do two matrix multiplies
     # both are ridiculous!
 
     #there's too many choices!!
     
     # if we focus on U, then we can either:
     
     # to get V:
     # -- compiute U.inv,  compute V = tcrossprod(U.inv)
     # -- compute W = crossprod(U), V = solve(W)
     # ^   IF we can find a way to do triangular-multiplies then the first should be cheaper by half a matrix; more numerically stable too
     # to get X:
     #  use backsolve
     
  # what. wait.. is this a PRECISION matrix or a COVARIANCE matrix?
  
  # "he Wishart distribution arises as the distribution of the sample covariance matrix" <https://en.wikipedia.org/wiki/Wishart_distribution>
  # "over the precision, \Lambda, is the Wishart distribution." <gaussian_prior_cheat_sheet.pdf>
  #  -> and this pdf very clearly says that the normal is sampled
  # but we SHOWED by experiment in HW1 that Wishart samples converge to the given COVARIANCE matrix, not its INVERSE
  #   --> is there some magic going on with the priors where you flip a precision matrix to a covariance matrix and vice versa?	
  # 
  
  # now, to gen X, we need the factorization of V:
  # U'U = W, which "sample" (actually we only sample A'A and then use hax)
  # V = W^-1
  # V = U^-1 * U'^-1
  # so U^-1 * z ~ N(0, U^-1 * U'^-1) = N(0, V)
  #  
  
   # on the other hand, if we decide we want to focus on U^-1
     #(assumption: no matter whether we use gamma or gamma.inv it isn't expensive because it isn't in the loop)
     # then we could either:
     #   -use U as above and invert
     # find A.inv and use
     #   - A.inv = backsolve(A, I)
     #   - U.inv = gamma.inv %*% A.inv #this... or:
     # then, to get V:
     #   - ?
     # to get X:
     #
  

BartlettFactor <- function(d, df) {
  # sample a cholesky factor of a standard Wishart distribution 
  #
  # args:
  #   d: the dimensionality
  #   df: the degrees of freedom of the Wishart distribution
  #
  # returns:
  #   a dxd matrix A ~ W(I_d, df)
  #
  # the Bartlett decomposition of M ~ W(I, df) is
  #   M = A'*A where 
  #   A is defined as an upper triangular matrix where
  #       A[i,i] ~ sqrt(X^2_{df-{i-1}})   
  #       A[i,j] ~ N(0,1)   (i<j)
  #
  #  it's not hard to derive from linearity of Normals and the def'n of the Wishart Dist that
  # CWC' ~ W(CIC', df), so you can get any Wishart you want from W(I, df)
  # (tho maybe it is hard to see with non-integer df. :shrug:)
  # so to extend this function to one which gives you W(V, df), do
  #    crossprod(chol(V)%*%BartlettFactor(d, df)) #XXX<-- double check this comment; it might be tcrossprod
  
  A = matrix(nrow=d, ncol=d)
  A[,] = 0
  df = df-(1:d)+1 #rchiqsq will vectorize over its DEGREES OF FREEDOM which is super clever
  diag(A) = sqrt(rchisq(d, df))
    
  # set the below-diagonals to N(0,1)s
  i_lower = col(A) > row(A) #XXX come out upper (EXPERIMENTAL)
  A[i_lower] = rnorm(sum(i_lower)) #sum() on a logical finds the len() of the lower triangle
  
  A
}



rNIW.version1 <- rNIW.typechecked(function(n, Mu, kappa, Psi, df) {
  # attempt to write out the efficient algorithm (but in R, of course)
#(one thing this setup doesn't do is precompute inv(Psi) ahead of time)
  d = length(Mu)
  ans = rNIW.alloc(n,d)

  # lawl wat: the canconical way to make an identity matrix in R is "diag(scalar)". wowwww.
  I = diag(d)

  #gamma.inv = chol(Psi) # upper triangular #<-- WRONG; I thought chol(inv(S)) = inv(chol(S)) but that's very not true
  gamma.inv = solve(chol(solve(Psi)))
  
  for(i in 1:n) {    
    # note: actually getting gamma = solve(gamma.inv) is expensive and numerically unstable,
    #       so you should avoid computing it if possible
    
    # construct the cholesky decomposition of a W(I, df) 
    A = BartlettFactor(d, df)
    #message("Upper triangular bartlett factor") #DEBUG
    #print(A)
    
    A.inv = backsolve(A, I)
    U.inv = gamma.inv %*% A.inv  #(U = AG and U'U = G'A'AG = W the Wishart-distributed matrix we never actually compute)
    
    ans$V[,,i] = tcrossprod(U.inv) #note well that this is not a LU decomposition -- it's UL
    
    # now we want X ~ N(Mu, V/kappa)
    z = rnorm(d); #sample N_d(0, I); since the covariance matrix is I, all draws are i.i.d. , so we can just sample d-univariates and reshape them into a vector
    ans$X[,i] = Mu + U.inv %*% z/sqrt(kappa)
  }
  
  ans
})


rNIW.version2 <- rNIW.typechecked(function(n, Mu, kappa, Psi, df) {
  # The difference between version1 and version2 is analogous to the difference between extremelynaive and naive:
  #  we factor out the common terms (Mu, kappa, chol(Psi)) to before and after the sampling loop
  
  d = length(Mu)
  ans = rNIW.alloc(n, d)
  
  #gamma.inv = chol(Psi)  #XXX Daniel proved this incorrect; leaving in now until can demonstrate that empirically with the testing code
  gamma.inv = solve(chol(solve(Psi)))
  
  I = diag(d)  # identity matrix; used for inverting via backsolve (solve(A) does inversion, but backsolve(A) says "missing argument"; sigh)
  for(i in 1:n) {  #apply scaling after the fact
     
     A = BartlettFactor(d, df)
     A.inv = backsolve(A, I)
     
     z = rnorm(d); #sample N_d(0, I); since the covariance matrix is I, all draws are i.i.d. , so we can just sample d-univariates and reshape them into a vector
     
     U.inv = gamma.inv %*% A.inv
     ans$V[,,i] = tcrossprod(U.inv)
     ans$X[,i] = U.inv %*% z
  }
  # in one step, scale everything
  # this relies on R's auto-splaying
  ans$X =  Mu + ans$X/sqrt(kappa) 
  
  ans
})

rNIW.version3 <- rNIW.typechecked(function(n, Mu, kappa, Psi, df) {
  #  instead of generating U.inv and using %*%, generate U and use backsolve
  # costs more at startup; perhaps that is much offset at runtime
  # XXX but because we need to output V, we necessarily need to construct U.inv anyway
  
  d = length(Mu)
  ans = rNIW.alloc(n, d)
  
  I = diag(d)  # identity matrix; used for inverting via backsolve (solve(A) does inversion, but backsolve(A) says "missing argument"; sigh)
  
  #gamma.inv = chol(Psi) #XXX daniel proved this wrong ; good thing, it's unused anyway
  #gamma = backsolve(gamma.inv, I)
  gamma = chol(solve(Psi))
  
  for(i in 1:n) {  #apply scaling after the fact
     
     A = BartlettFactor(d, df)
     
     U = A %*% gamma
     
     # sample V
     U.inv = backsolve(U, I) 
     ans$V[,,i] = tcrossprod(U.inv)
     
     # sample X
     z = rnorm(d); #sample N_d(0, I); since the covariance matrix is I, all draws are i.i.d. , so we can just sample d-univariates and reshape them into a vector
     ans$X[,i] = backsolve(U, z)
  }
  # in one step, scale everything
  # this relies on R's auto-splaying
  ans$X =  Mu + ans$X/sqrt(kappa) 
  
  ans
})



# TODO: try replacing every %*% with triangularmultiply (which is available in LAPACK as dtrmm)
#            backsolve a form of triangularmultiply, which does two operations at once:
#              solve(T)%*%v == backsolve(T, v)

# TODO: try using BLAS routines; the canonical interface to this in R is the "Matrix" package
#   and that package has all sorts of shades:
# TODO: try using all different types of triangular matrices and seeing how the speed changes:
#  dtrMatrix - "dense" format (all values stored, but the matrix knows its a triangular matrix and ops will only operate on half of it)
#  dtRMatrix - "sparse" format
#  dtCMatrix - "compressed sparse" format
#  dtpmatrix - "packed" format (only the d(d+1)/2 non-null values are stored in one long vector and indexing tricks are used to make this okay
#TODO: write a TriBartlettFactor() which returns a dtpMatrix
# ^ THE ABOVE IS ALL FULL OF FAIL: the Matrix package uses the S4 object system
 #, and S4's lookup routines have huge overhead
# It might be worth it when d and n grow very very large
# and it would even be interesting to code them up and see where that point is
# but for now the idea is scrapped


#TODO: figure out if there is a clever linear algebra trick which would let us
#            multiply gamma across every sample in one step
#            this will not give an improvement on standard CPU, but if we had a GPU handy it might be a giant boon.


# hand-rolled C algorithm  
rMNIW.Rcpp <- rMNIW.typechecked(function(n, Mu, Kappa, Psi, df){
  
  # precompute some values that would be a headache to do in C
  d = dim(Mu)[2];
  q = dim(Mu)[1];
  
  Kappa_Inv <- chol(solve(Kappa)); # to give sqrt of kappa_Inv later
  gamma.inv = solve(chol(solve(Psi)));
  return(rMNIW_Rcpp(n,d,q,Mu,Kappa_Inv, gamma.inv, df));
})


rNIW.Rcpp <- rNIW.typechecked(function(n, Mu, kappa, Psi, df) {
  # precompute what can be precomputed
  # because the slowness of doing them in R will not be
  # that much (only O(1), not O(n)), and is outweighed by the headache in C
  # TODO: do these in C as well, for completeness.
  
  # precompute some useful values
  d = length(Mu)
  gamma.inv = solve(chol(solve(Psi)))
    
  return(rNIW_Rcpp(n, d, Mu, kappa, gamma.inv, df))
})


# efficient algorithm using the Eigen library
rNIW.version.RcppEigen <- rNIW.typechecked(rNIW_version_eigen)


rNIW.based_on_multi <- rNIW.typechecked(function(n, Mu, Kappa, Psi, df) { # TODO: pick a better name
  # rNIW implementation backed by rMNIW
  # TODO: fill this in and kill rNIW_Rcpp
  # you can view this as a convenience wrapper
  # 
  
  # coerce Mu and Kappa to matrices 
  d = length(Mu)
  dim(Mu) = c(1, d)
  dim(Kappa) = c(1,1)
  
  rMNIW(n, Mu, Kappa, Psi, df)
})

# Final decision about the exported implementation:
rNIW = rNIW.Rcpp
rMNIW = rMNIW.Rcpp
