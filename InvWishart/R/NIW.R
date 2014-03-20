
require(MASS); #for rmvnorm
rmvnorm <- mvrnorm; mvrnorm <- NULL; #repair naming convention

require("Matrix") #for faster(?) solving routines (this package is very slow to load; comment this line out if you must)
                  #WARNING: with this loaded the S4-methodlookup overhead seems to drastically slow down rNIW.extremelynaive, rNIW.naive and rNIW.snappy1
##################################
### Naive Sampling

rInvWishart <- function(n, df, Sigma) {
 # inverse wishart sampler 
 stopifnot(nrow(Sigma) == ncol(Sigma))
 stopifnot(df >= nrow(Sigma))
 W = rWishart(n, df, solve(Sigma));
 for(i in 1:n) {
   W[,,i] = solve(W[,,i]) 
 }
 W
}


rNIW.extremelynaive <- function(n, d, Mu, kappa, Psi, df, V, X) {
  # generate the n samples by n times doing: i) generate V, ii) generate X given V.

  for(i in 1:n) {
     V[,,i] = rInvWishart(1, df, Psi)[,,1]; #sample only one 
     X[,i] = rmvnorm(1, Mu, V[,,i]/kappa);  #ditto
  }
  list(V=V, X=X)
}


rNIW.naive <- function(n, d, Mu, kappa, Psi, df, V, X) {
  # generate the n samples
  # slightly less totally durpy implementation
  # the main difference between this and extremelynaive is that this only inverts Psi once
  
  V = rInvWishart(n, df, Psi);   #ignore the prealloc'd space
  
  # must be in a loop because each sample has a different distribution
  for(i in 1:n) {
     X[,i] = rmvnorm(1, Mu, V[,,i]/kappa);
  }
  list(V=V, X=X)
}




#######################################################
### SnappySampling by exploiting Bartlett Decomposition

BartlettFactor <- function(d, df) {
  # the Bartlett decomposition of X ~ W(V, df) is
  #   X = gamma*A*A'*gamma'
  # where  L = chol(V) (lower-triangular chol) and
  #        A is defined as A[i,i] ~ sqrt(X^2_{df-{i-1}}) && A[i,j] ~ N(0,1) and A[j,i] = 0 (i>j)
  
  # the base case (or a special case, depending on how you look at it) is
  #  W = AA' which has W ~ W(I, df)  <-- this is the Bartlett Theorem
  #  it's not hard to derive from linearity of Normals and the def'n of the Wishart Dist that
  # CWC' ~ W(CIC', df), so you can get any Wishart you want from W(I, df)
  # (tho maybe it is hard to see that with non-integer df. :shrug:)
  
  A = matrix(nrow=d, ncol=d)
  A[,] = 0
  df = df-(1:d)+1 #rchiqsq will vectorize over its DEGREES OF FREEDOM which is super clever
  diag(A) = sqrt(rchisq(d, df))
    
  # set the below-diagonals to N(0,1)s
  i_lower = col(A) > row(A) #XXX come out upper (EXPERIMENTAL)
  A[i_lower] = rnorm(sum(i_lower)) #sum() on a logical finds the len() of the lower triangle
  
  A
}


rNIW.snappy1.sample <- function(d, Mu, kappa, Psi, df) {
  # attempt to write out the efficient algorithm (but in R, of course)
  
  
  gamma.inv = chol(Psi) # upper triangular
  # note: actually getting gamma = solve(gamma.inv) is expensive and numerically unstable,
  #       so you should avoid computing it if possible
  
  # construct 
  A = BartlettFactor(d, df)
  #message("Upper triangular bartlett factor")
  #print(A)
  
  # okay, well, it seems we're stuck with inverting A
  U = gamma
  # lawl wat: identity matrix in R is "diag(scalar)". wowwww.
  A.inv = backsolve(A, diag(d))
  U.inv = gamma.inv %*% A.inv  #(U = AG and U'U = G'A'AG = W the Wishart-distributed matrix we never actually compute)
  
  V = U.inv %*% t(U.inv) #note well that this is not a LU decomposition -- it's UL
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
  
  # now we want X to follow some stuff
  # one method (slow and numerically unstable)
  z = rnorm(d); #sample N_d(0, I); since the covariance matrix is I, all draws are i.i.d. , so we can just sample d-univariates and reshape them into a vector
  X = Mu + U.inv %*% z/sqrt(kappa)
  
list(X=X, V=V);
}

#(one thing this setup doesn't do is precompute inv(Psi) ahead of time)

# 



rNIW.snappy1 <- function(n, d, Mu, kappa, Psi, df, V, X) {
  # generate the n samples by n times doing: i) generate V, ii) generate X given V.
  for(i in 1:n) {
     sample = rNIW.snappy1.sample(d, Mu, kappa, Psi, df)
     V[,,i] = sample$V
     X[,i] = sample$X
  }
  list(V=V, X=X)
}





rNIW.snappy2.sample <- function(d, df) {
  A = BartlettFactor(d, df)
  A.inv = backsolve(A, diag(d))
    
  z = rnorm(d); #sample N_d(0, I); since the covariance matrix is I, all draws are i.i.d. , so we can just sample d-univariates and reshape them into a vector
  X = A.inv %*% z
list(X=X, A.inv=A.inv);
}

rNIW.snappy2 <- function(n, d, Mu, kappa, Psi, df, V, X) {
  # generate the n samples by n times doing: i) generate V, ii) generate X given V.
  # analogous to the difference between extremelynaive and naive,
  # the only difference from snappy1 is that we only factor Psi once
  
  gamma.inv = chol(Psi) 
  kappy = sqrt(kappa)
  for(i in 1:n) {  #apply scaling after the fact
     sample = rNIW.snappy2.sample(d, df)
     U.inv = gamma.inv %*% sample$A.inv   #this seems dumb. TWO multiplies by gamma? hm. if 
     V[,,i] = U.inv %*% t(U.inv)
     X[,i] = Mu + gamma.inv%*%sample$X/kappy 
  }
  list(V=V, X=X)
}


# is it better to actually construct U than try to apply it after the fact?
#   --> if all we wanted was X, then no, but we also need V, which means we MUST construct U somewhere
# (notice: U is the factor that is necessary...)

# snappy3: apply Mu, gamma and kappa after the generating loop
#     hope: vectorization will kick in and make this superfast


rNIW.snappy3 <- function(n, d, Mu, kappa, Psi, df, V, X) {
  # generate the n samples by n times doing: i) generate V, ii) generate X given V.
  # analogous to the difference between extremelynaive and naive,
  # the only difference from snappy1 is that we only factor Psi once
  gamma.inv = chol(Psi) 
 
  for(i in 1:n) {  #apply scaling after the fact
     sample = rNIW.snappy2.sample(d, df) #USING SNAPPY2 STILL
     U.inv = gamma.inv %*% sample$A.inv
     V[,,i] = U.inv %*% t(U.inv)
     X[,i] = sample$X
  }
  
  # now, how do I apply a matrix multiply to a whole vector at once?
  # ..I should just be able to say
  X = Mu + gamma.inv %*% X / sqrt(kappa) # because X has samples column-wise in this setup (the index i is in the last component)
  
  # is there a way to do the same for V?
  # probably not... 
  
  list(V=V, X=X)
}


# snappy 2.5 replaces %*% with appropriate slightly-more-efficient matrix multiplies:
#      a "symmetrictriangularmultiply" (aka dsyrk.f, but if we can avoid using fortran that would be good)
#      a "upper triangle times upper triangle" multiply
# this requires using the Matrix package

rNIW.snappy2.5.sample <- function(d, df) {
  A = BartlettFactor(d, df)
  A= as(A, "dtrMatrix") #cast so we can invoke LAPACK routines
  A.inv = solve(A)  #Matrix knows how to invert itself by backsolving
                    #help(Matrix::solve) --> '          all use LAPACK based versions of efficient triangular
                    #                        'backsolve', or 'forwardsolve'.'
  
  z = rnorm(d); #sample N_d(0, I); since the covariance matrix is I, all draws are i.i.d. , so we can just sample d-univariates and reshape them into a vector
  X = A.inv %*% z
list(X=X, A.inv=A.inv);
}


rNIW.snappy2.5 <- function(n, d, Mu, kappa, Psi, df, V, X) {
  Psi = Matrix(Psi); #cast to Matrix package matrices, which provide C/Fortran case-specific efficient matrix ops
  Mu = Matrix(Mu);

  gamma.inv = chol(Psi)
  kappy = sqrt(kappa)
  for(i in 1:n) {  #apply scaling after the fact
     message(i)
     sample = rNIW.snappy2.5.sample(d, df)
     U.inv = gamma.inv %*% sample$A.inv   #this seems dumb. TWO multiplies by gamma? hm. 
     V[,,i] = matrix(tcrossprod(U.inv)) #compute U.inv %*% t(U.inv)  
     X[,i] = matrix(Mu + gamma.inv%*%sample$X/kappy )  #are the matrix() calls slowing us down??
  }
  list(V=V, X=X)
}

# TODO: only use U instead of using A in one place and U in another

# TODO: try using all different types of triangular matrices and seeing how the speed changes:
#  dtrMatrix - "dense" format (all values stored, but the matrix knows its a triangular matrix and ops will only operate on half of it)
#  dtRMatrix - "sparse" format
#  dtCMatrix - "compressed sparse" format
#  dtpmatrix - "packed" format (only the d(d+1)/2 non-null values are stored in one long vector and indexing tricks are used to make this okay

#TODO: try making the Bartlett factor matrix triangular from the get-go (we could even use a dtpMatrix..)

# snappy 2.7 additionally uses a better backsolve

# the diff with snappy4 is that it does
# backsolve() + forwardtriangularmultiply(???) instead of computing U.inv


## snappy5 has figured out how to multiply gamma across things in one step

## snappy6 eschews BartlettFactor and instead runs InverseBartlettFactor
## (getting to this point will require some very very careful derivations that may be beyond us)

#########################3
## Hosting code

rNIW.typechecks <- function(rNIW.sampleloop, n, Mu, kappa, Psi, df) {
  # n: the desired number of sample points (X, V) ~ NIW(Mu, kappa, Psi, df)
  # df need not be integer; instead it can be anything in the domain of the gamma function
  #  rNIW.sampleloop: a continuation (as in continuation-passing-style) to go to when the typechecks are done
  
  # this subroutine is just the safe(ish) type-checked outer call
  # which calls rNIW.sample to get the real samples (so you can experiment with swapping out different samplers)
  
  Mu = as.vector(Mu)
  Psi = as.matrix(Psi)
  
  d = dim(Psi)[1]
  stopifnot(d == dim(Psi)[0]) #Psi must be square
  stopifnot(d == length(as.vector(Mu)))  #Mu must be the same dimensionality as the covariance matrix it supposedly goes with

  
  # pre-allocate space
  X = rep(0, d*n);
  dim(X) = c(d, n);

  V = rep(0, d^2*n);
  dim(V) = c(d,d,n);

    
  rNIW.sampleloop(n, d, Mu, kappa, Psi, df, V, X)
}


################################3
### TESTING


kMu = c(0.25, -1.5, 0.33, 9)
kKappa = 1
kV  = cbind(c(2.84, 0.43, 0.16, .5), c(0.43, 1.52, -0.24, -.73), c(.16, -.24, 4.49, 88), c(1,-2,3,4))
kV  = t(kV)%*%kV #after I added digits at random, force kV to be some positive def matrix
kDF = 7.32 # must be at least as large as the dimension of kV


test <- function(method, N=32324) {

  tryCatch({ #wrapped in a tryCatch so that tests can fail and the code can plug along anyway

  name = deparse(substitute(method)) #LOL R


message(paste("Beginning",name,"..."))
tic = proc.time()
Samples = rNIW.typechecks(method, N, kMu, kKappa, kV, kDF)
toc = proc.time()
runtime = toc - tic

# time it!
message(paste("Runtime for", name,":"))
print(runtime)


# plot it. plot it alllll.

par(mfrow=c(1,1))
plot(c(), c(), xlim=c(0,0), ylim=c(0,0), main=name, sub="you love it")

par(mfrow=c(2,2))
for(i in 1:4) { for(j in 1:4) { #XXX hardcoded
  hist(Samples$V[i,j,], main=paste("Hist of V[",i,",",j,"]"),  probability=T, breaks=20, col='red')
  if(exists("Baseline.Samples")) {
    lines(density(Baseline.Samples$V[i,j,]), col='black')
  }
}}

for(i in 1:4) {
  hist(Samples$X[i,], probability=T, breaks=20, col='red', main=paste("Hist of X[",i,"]"))
  if(exists("Baseline.Samples")) {
    lines(density(Baseline.Samples$X[i,]), col='black')
  }
}

message()
Samples #return for future use, in case you care
}, 
 error = function(e) { print(e); },
 warning = function(e) { print(w); }
 )
}  


message()
message("Starting test runs")
message("------------------")
#Baseline.Samples = test(rNIW.extremelynaive)
#ignored = test(rNIW.naive)
#ignored = test(rNIW.snappy1);
#ignored = test(rNIW.snappy2);
ignored = test(rNIW.snappy2.5);
#ignored = test(rNIW.snappy3);
