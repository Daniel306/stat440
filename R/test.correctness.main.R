# test.correctness.main.R
#
#  Until we learn the proper Rish way of writing tests.
#  This script will serve.
#  in Python, the contents of this file would be under the
#   "if __name__ == '__main__'" 
#  block.

source("test.correctness.R")


################################3
### TESTING


kMu = c(0.25, -1.5, 0.33, 9)
kKappa = 1
kPsi  = cbind(c(2.84, 0.43, 0.16, .5), c(0.43, 1.52, -0.24, -.73), c(.16, -.24, 4.49, 88), c(1,-2,3,4))
kPsi  = t(kV)%*%kV #after I added digits at random, force kV to be some positive def matrix
kDF = 7.32 # must be at least as large as the dimension of kV


test <- function(method, N=32343) {
  # XXX THIS IS BROKEN AFTER REFACTOR
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


#rNIW.naive = wrap(rNIW.naivecore) somethiung something

ground = rNIW.naive(500, kMu, kKappa, kPsi, kDF)

samples = (function(n, mu, kappa, psi, df) {
  d = length(mu)
  X = matrix(NA, d, n)
  V = matrix(NA, d*d, n)
  dim(V) = c(d,d,n)
  rNIW.snappy3.5(n, d, mu, kappa, psi, df, X, V)
})
ignored = test(rNIW.snappy3.5);


# our tests of interest:
# 1) sample values
#  -> marginal densities
#  -> moments
# 2) runtimes
# these purposes feel like they should share code
# but do not be overzealous: 
# runtimes 
