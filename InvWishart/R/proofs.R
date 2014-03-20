
# test various properties that may/may not help us

n = 5
d = 3

X = rnorm(n*d)
dim(X) = c(n, d)
#print(X) #a totally non-symmetric matrix

message("S")
S = t(X)%*%X #guaranteed to be symmetric
print(S)

message("L")
L = t(chol(S)) #NB: R's chol() gives an upper triangle instead of lower
print(L)

# sanity check: LL' = S
stopifnot(norm(L%*%t(L) - S) < 1e-10)

# 

#------------


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

R = sapply(1:222, function(i) { solve(BartlettFactor(5, 6)) })
# sapply flattens the results into a single vector
# but we can put it back like this:
dim(R) = c(sqrt(dim(R)[1]), sqrt(dim(R)[1]), dim(R)[2])

# and now we can check out the univariate dist of any single BartlettFactor entry like..
x = seq(-30, 30, 0.01)
hist(R[1,1,], probability=T, breaks=20)
lines(x, dchisq(x, 6), col='blue')

hist(R[3,3,], probability=T, breaks=20)
lines(x, dchisq(x, 6-3+1), col='blue')


hist(R[2,5,], probability=T, breaks=200)
lines(x, dnorm(x), col='blue')


R.cor = rep(NA, dim(R)[1]^4)
dim(R.cor) = c(dim(R)[1], dim(R)[1], dim(R)[1], dim(R)[1])
for( i in 1:dim(R)[1]) { for( j in 1:dim(R)[2]) {
for(i2 in 1:dim(R)[1]) { for(j2 in 1:dim(R)[2]) {
  R.cor[i,j,i2,j2] = cor(R[i,j,], R[i2,j2,])
}}
}}

print(mean(na.omit(as.vector(R.cor > .10))))


# okay, now what happens if we INVERT
for(k in 1:dim(R)[3]) {
  R[,,k] = solve(R[,,k])
}

# look for correlations by brute force
# even though the inverse involves all columns,
# maybe the entries are (effectively?) uncorrelated
# and it is just as good to compute them directly
# e.g. using F distributions and stuff
R.cor = rep(NA, dim(R)[1]^4)
dim(R.cor) = c(dim(R)[1], dim(R)[1], dim(R)[1], dim(R)[1])
for( i in 1:dim(R)[1]) { for( j in 1:dim(R)[2]) {
for(i2 in 1:dim(R)[1]) { for(j2 in 1:dim(R)[2]) {
  R.cor[i,j,i2,j2] = cor(R[i,j,], R[i2,j2,])
}}
}}

print(mean(na.omit(as.vector(R.cor > .10))))

# ^ huh. 
# the inverse is consistently LESS (linear-pearson) correlated than the original
