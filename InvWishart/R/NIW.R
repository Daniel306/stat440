
require(MASS); #for rmvnorm
rmvnorm <- mvrnorm; mvrnorm <- NULL; #repair naming convention


rInvWishart <- function(n, df, Sigma) {
 # inverse wishart sampler 
 W = rWishart(n, df, solve(Sigma));
 for(i in 1:n) {
   W[,,i] = solve(W[,,i]) 
 }
 W
}

rNIW.naivesample <- function(Mu, kappa, Psi, df) {
  # naive slow-as-mud Normal-Inverse-Wishart sampler
  V = rInvWishart(1, df, Psi)[,,1]; #sample only one 
  X = rmvnorm(1, Mu, V/kappa);
  
  list(X=X, V=V);
}


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
  i_lower = col(A) < row(A)
  A[i_lower] = rnorm(sum(i_lower)) #sum() on a logical finds the len() of the lower triangle
  
  A
}

rNIW.snappysample <- function(Mu, kappa, Psi, df) {
  d = length(Mu) #just shorthand, not for typechecking
  # attempt to write out the efficient algorithm (but in R, of course)
  
  
  gamma = t(chol(Psi)) #t() because in R chol gives the upper triangle (the right term) and I wrote my derivations in terms of the lower
  
  # 
  
     # construct L_i
  L_i = BartlettFactor(d, df)
  message("Lower triangular bartlett factor")
  print(L_i)
  
  # so we have that LL' ~ W(I, df) (or, alternately, U'U = LL')
  # we *want* V ~ W^-1(Psi, df)
  #          <=> inv(V) ~ W(inv(Psi), df)
  #              inv(V) ~ W(chol(inv(Psi))' * chol(inv(Psi))), df)
  #         if we define inv(V) = [ chol(inv(Psi))' * U' ] * [ U * chol(inv(Psi)) ] then inv(V) ~ W(inv(Psi), df)
  #                  so  inv(inv(V)) ~ inv(W)(Psi, df)
  #           hurrah!
  #         now how to actually ~use~ this...
  # 
  #
  # we have gamma = chol(Psi)'; #lower triangular
  #       chol(inv(Psi))' = inv(chol(Psi)) = inv(gamma') = inv(gamma)'
  #  and  chol(inv(Psi))  = inv(chol(Psi))'= inv(gamma')'= inv(gamma)'' = inv(gamma)
  #   so inv(V) = inv(gamma)' * U' * U * inv(gamma)
  #  and thus V = inv(inv(gamma)' * U' * U * inv(gamma)) = {note the reversal in order} gamma * inv(U) * inv(U)' * gamma'
  # i'm making some really simple mistake here.. maybe gamma isn't supposed to be lower on the left?
  # yes! that's it! gamma is supposed to be UPPER on the left
  
  U_i = gamma%*%t(solve(L_i)) #hmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
  print(U_i)
  si
  
  # 
  # now, Sb%*%t(Sb) ~ W(I, n) according to Bartlett's Method according to the assignment
  # but we want W(V, n-1)
  # we can use part i), which says if S ~ W(I, n) then LSL' ~ W(LIL', n)  (NB, in case i forget: ' means "transpose")
  # We want LIL' = LL' to be V, so we want to apply the choklesky decomposition to get square roots:
  L = chol(V) # note: R's chol() gives an *upper* triangle... 
    # and L'L = V; tho since the inner term of W(I, n) is I, it doesn't matter if we decompose the wrong way..
    # since S ~ W(I, n) => L'SL ~ W(L'IL, n) ~ W(V, n)
  # Sb%*%t(Sb) ~ W(I, n), so (L%*%Sb)%*%t(L%*%Sb) = L(SbSb')L' ~ W(LL', n)
  # so we can build Sb subsequently by taking the lower triangular matrix, scaling, and squaring:
  Sb = t(L)%*%Sb  #lower-triangular * lower triangular = lower-triangular
  Sb = Sb%*%t(Sb)
  Sb

  
     # construct W_i (?)
     # construct U_i (?)
     # construct V_i

list(X=X, V=V);
}

#(one thing this setup doesn't do is precompute inv(Psi) ahead of time)

# 

rNIW.sample = rNIW.naivesample
rNIW.sample = rNIW.snappysample

rNIW.wrapper <- function(n, Mu, kappa, Psi, df) {
  # n: the desired number of sample points (X, V) ~ NIW(Mu, kappa, Psi, df)
  # df need not be integer; instead it can be anything in the domain of the gamma function
  
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
  
  # generate the n samples
  for(i in 1:n) { 
     message(paste("sample #", i, sep=""))
     sample = rNIW.sample(Mu, kappa, Psi, df);
     X[,i] = sample$X
     V[,,i] = sample$V
     message("X:")
     print(sample$X)
     message("V:")
     print(sample$V)
     message()
  }
  list(V=V, X=X)
}


rNIW = rNIW.wrapper

kMu = c(0.25, -1.5, 0.33, 9)
kKappa = 1
kV  = cbind(c(2.84, 0.43, 0.16, .5), c(0.43, 1.52, -0.24, -.73), c(.16, -.24, 4.49, 88), c(1,-2,3,4))
kV  = t(kV)%*%kV
kDF = 3.3


Z = rNIW(6, kMu, kKappa, kV, kDF)

#print(Z)
