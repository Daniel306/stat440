

# convention: 'X.typecheck' has the preconditions (written as assertions, which will crash if the caller gives wrong data) that actually do validation.
#                   'X.typechecked' wraps a function implementing API "X" with "X.typecheck" and additionally checks postconditions.

NIW.typecheck <- function(Mu, Kappa, Psi, df) {
    # Common-core of typechecks of the NIW distribution's parameters
    #
    # These are things like "Psi is symmetric" and "Mu has the same dimensionality as Psi".
    # TODO: this clearly shares DNA with MNIW.typecheck, but; try harder to factor it
    
    # the Xs in (X,V)~NIW(...) are "univariate": each X is one column:
    stopifnot(is.vector(Mu))
    d = length(Mu)
    
    # as a consequence, Kappa must be scalar
    stopifnot(length(Kappa) == 1)
    
    # Kappa must be positive definite
    # because Kappa is scalar for NIW, we can actually check this directly:
    stopifnot(Kappa > 0)
    
    #Psi must be positive definite
    # XXX we don't actually check that Psi is positive definite because that's hard and slow
    #  we check for symmetry (and squareness) which is relatively cheap 
    stopifnot(is.symmetric(Psi)) #implies "is.matrix()"

    #Mu must be the same dimensionality as the covariance matrix it supposedly goes with
    stopifnot(nrow(Psi) == length(Mu))
    
    # df is a scalar > d-1
    stopifnot(length(df) == 1)
    stopifnot(df > d - 1)
}

MNIW.typecheck <- function(Mu, Kappa, Psi, df) {
    # Matrix-Normal|Inverse-Wishart typechecks
    # setup / motivation
    # in Y = XB + E ~ MN(Mu, I, Kappa)
    # let q be the number of response variates (columns of Y)
    # let d be the number of predictor variates (columns of X)
    # then
    #  1) Mu should be a q x d matrix
    #  2) Kappa should be a q x q matrix (this gives the covariance between _____________?)
    #  3) Psi should be a d x d matrix (this gives **the prior** for the covariance between __________?)
    #  4) Kappa and Psi should be positive definite (XXX this isn't actually checked)
    
    # the Xs in (X,V)~NIW(...) are multivariate: each X is a whole matrix
    stopifnot(is.matrix(Mu))
    q = dim(Mu)[1]
    d = dim(Mu)[2]
    
    # Kappa must be positive definite
    stopifnot(is.symmetric(Kappa))
    # and the same size as
    stopifnot(dim(Kappa)[1] == q)
    
    #Psi must be positive definite
    # XXX we don't actually check that Psi is positive definite because that's hard and slow
    #  we check for symmetry (and squareness) which is relatively cheap 
    stopifnot(is.symmetric(Psi)) #implies "is.matrix()"
    #Mu must be the same dimensionality as the covariance matrix it supposedly goes with
    stopifnot(dim(Psi)[1] == d)
    
    # df is a scalar > d-1
    stopifnot(length(df) == 1)
    stopifnot(df > d - 1)
}
