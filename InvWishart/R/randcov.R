
#

# 


randcov = function(d) {
# generate random positive semidefinite matrices
#  -
# http://stats.stackexchange.com/questions/2746/how-to-efficiently-generate-positive-semi-definite-correlation-matrices
# - http://mail.scipy.org/pipermail/scipy-user/2007-June/012487.html
# - http://stackoverflow.com/questions/1040324/how-to-generate-pseudo-random-positive-definite-matrix-with-constraints-on-the-o

# this trick is probably good enough for our purposes:
  rWishart(1, d+1, diag(d))[,,1]
}


# time how long the cholesky decomposition takes
timechol <- function() {
  R = data.frame(d=c(), time=c())
  for(d in 1:3200) {
    message("d=",d)
    for(sample in 1:35) {
      Psi = randcov(d);
      #message("time it")
      z = proc.time()
      U = chol(Psi)
      z = proc.time() - z
      #message("record it")
      R = rbind(R, c(d, z["elapsed"]))       #LAWL SLOW
    }
    write.csv(R, "randcov.times.Rdata")
  }
  R
}

R = timechol()

