
# timings on R's BLAS interface "Matrix"




###################3
### instrumentation

timeit <- function(thunk, n) {
  times = rep(-171717, n);
  message("hai");

  for(i in 1:n) {
    message(i);
    tryCatch( {
      ptm = proc.time();
      thunk();
      times[i] = (proc.time() - ptm)["elapsed"]; #record elapsed time
    }, error = function(e) { print(e) }) #TODO: log errors
  }
  times;
} 

plotit <- function(thunk, n) {
  title = paste(deparse(substitute(thunk)))
  measurement = timeit(thunk, n)
  message(title)
  print(measurement)
  hist(measurement, main=title)
  plot(density(measurement), main=title)
  #plot(density(measurement, kernel="rectangular"), main=title)
  message()
}


#####################3 test

#require("Matrix")

#n = 2411; #this number doesn't actually matter that much except for stabilizing chol()
d = 333;

# bah, this code is clever but does not work because I don't have a way
# to generate invertible upper triangular matrices.

b = rnorm(d); #target vector: the test task is solving y such that Uy = b (ie compute inv(U)*b)

make_upper_triangle <- function() {
  # alg3: <http://www.ams.org/notices/200705/fea-mezzadri-web.pdf>
  # not totally implemented as you can see; I didn't get the parts that stabilize to work, bc built in qr() in R is difficult
  X = rnorm(d^2)
  dim(X) = c(d,d)
  #print(X)
  X = (qr(X)$qr)
  #Q = X
  R = X
  R[row(X) > col(X)] = 0 #kill lower triangle of R
  #d = diag(R)
  #ph = d / abs(d)
  #print(Q)
  #print(R)
  #print(ph)
  #Sys.sleep(1)
  #X = Q %*% t(ph) %*% X
  R

  # alg2: (unstabllle)
  #X = matrix(ncol=d, nrow=d);
  #X[,] = 0
  #idx = row(X) <= col(X)
  #X[idx] = rnorm(sum(idx))
  #X
  
  # alg1: (slowwww)
  #stopifnot(n >= d);
  #X = rnorm(n*d);
  #dim(X) = c(n, d);
  #S = t(X)%*%X;
  #U = chol(S);
  #U
}

triangular.regular.regular <- function() {
  U = make_upper_triangle();
  solve(U, b)
}

triangular.regular.backsolve <- function() {
  U = make_upper_triangle();
  backsolve(U, b)
}

triangular.regular.blas <- function() {

}

triangular.sparse.blas <- function() {

}

n = 250
#plotit(make_upper_triangle, n)
plotit(triangular.regular.regular, n)
plotit(triangular.regular.backsolve, n)

warnings()
