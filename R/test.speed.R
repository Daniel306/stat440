
# something like

# TODO: NIW.runtime which only runs one algorithm at a time??
# and then .runtimes() which cbinds the results from that all together? maybe?
# hm

NIW.runtimes <- function(algorithms, N=10^(1:10) log10, rep=5) {
  # Record the runtimes of a suite of Normal-Inverse-Wishart samplers, for comparison.
  #
  # args:
  #  algorithms: a character vector containing the names of functions to call (which will be looked up with 'get()')
  #   each function must have this API:
  #     f(N)
  #   so to
  #  N:   the points at which to sample the runtimes algorithm (n being the number of samples to ask each NIW to generate)
  #  rep: the number of times to sample each runtime
  #
  # returns:
  #   the results in a matrix (n, algorithm, time)
  #   there should be 'rep' many rows for each pair (n, algorithm)
  #
  #  TODO:
  #   - also sample the runtimes across different values of the distribution paramters
  #       ((we expect that most runtimes should be constant in Mu, Kappa, V, and slightly increasing with DF
  #         and the naive algorithm should show bad degredation as length(Mu) == dim(V)[1] == dim(V)[2] grows))
 
  # (algorithm is a factor column)
  R = data.frame(matrix(NA, length(N)*length(algorithms)*rep, 3))
  colnames(R) = c("n","algorithm","time")

  # i is the index into R; really I wish for an "enumerate()" to loop over a crossproduct
  # but this will do
  i = 1;
  for(a in algorithms) { # These loops
    A = get(a)           # are flat
  for(n in N) {          # on purpose
  for(r in 1:rep) {      # to match the flat matrix output they are generating
    tic = proc.time()
    .unused.samples = A(n) # TODO: pass in kMu, kKappa, kV, kDF)
    toc = proc.time()

    R[i,] = c(n, a, (toc - tic))
    i = i + 1;
  }
  }
  }

  #coerce the algorithm column to the categorical variable it is
  # this needs to happen here because 'factors' are just thinly disguised enums
  # and there's no way to e.g. assign a string to a pre-existing factor
  R$algorithm = as.factor(R$algorithm)
  return(R)
}