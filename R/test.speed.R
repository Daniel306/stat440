

source("NIW.R")
source("test.constants.R")

# TODO: there's some appeal in a NIW.runtime() which only runs one algorithm at a time
#       and then writing NIW.runtimes() to just rbinds the results from that all together?

NIW.runtimes <- function(algorithms, max=1e7, N=floor(exp(seq(0, log(max), length.out=15))), rep=5) {
  # Record the runtimes of a suite of Normal-Inverse-Wishart samplers, for comparison.
  #
  # args:
  #  algorithms: a character vector containing the names of functions to call (which will be looked up with 'get()')
  #   each function has the rNIW API as documented in NIW.R
  #  N:   the points at which to sample the runtimes algorithm (n being the number of samples to ask each NIW to generate)
  #  max: the largest n to test, if N is left at its default (and unused otherwise)
  #  rep: the number of times to sample each runtime
  #
  # returns:
  #   the results in a matrix (n, algorithm, time)
  #   there should be 'rep' many rows for each pair (n, algorithm)
  #   and algorithm is a factor column.
  #
  #  TODO:
  #   - also sample the runtimes across different values of the distribution paramters
  #       ((we expect that most runtimes should be constant in Mu, Kappa, V, and slightly increasing with DF
  #         and the naive algorithm should show bad degredation as length(Mu) == dim(V)[1] == dim(V)[2] grows))
  #     currently just uses the values in test.constants.R
 
  R = data.frame(matrix(NA, length(N)*length(algorithms)*rep, 3))
  colnames(R) = c("n","algorithm","time")

  # i is the index into R; really I wish for an "enumerate()" to loop over a crossproduct
  # but this will do for now.
  i = 1;
  for(a in algorithms) {               # These loops are
    A = get(paste("rNIW", a, sep=".")) # purposely flat
  for(n in N) {                        # to match the flat matrix
  for(r in 1:rep) {                    # they generate.
    message("runtime a=",a,", n=", n, ", rep=",r)  #DEBUG
    tic = proc.time()
    A(n, kMu, kKappa, kPsi, kDF) #generate samples
    toc = proc.time()
    #print(toc - tic) #DEBUG

    # because R is a jerk
    # these lines need to be done one at a time
    # otherwise the presence of the string a makes the whole tuple type "character"
    # and that causes the whole dataframe to be "character" typed
    R[i,"n"] = n
    R[i, "algorithm"] = a
    R[i,"time"] = (toc - tic)["elapsed"]
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