# test.correctness.main.R
#
#  Until we learn the proper Rish way of writing tests.
#  This script will serve.
#  in Python, the contents of this file would be under the
#   "if __name__ == '__main__'" 
#  block.

source("test.correctness.R")

main <- function() { # avoid polluting the namespace

    # Distribution parameters for testing with
    # TODO: let these be user set and/or randomize them
    kMu = c(0.25, -1.5, 0.33, 9)
    kKappa = 1
    kPsi  = cbind(c(2.84, 0.43, 0.16, .5), c(0.43, 1.52, -0.24, -.73), c(.16, -.24, 4.49, 88), c(1,-2,3,4))
    kPsi  = t(kPsi)%*%kPsi #after I added digits at random, force kV to be some positive def matrix in the most crude way I know
    kDF = 7.32 # must be at least as large as the dimension of kV
    
    # number of samples to trial
    # note that ground uses a different number, because rNIW.naive is slow
    n = 32343
    
    message()
    message("------------------------------------")
    message("Generating reference samples.")
    ground = rNIW.naive(15500, kMu, kKappa, kPsi, kDF)
    #plot.NIW.marginals(ground, ground, alg="naive") #DEBUG
    #ks.test.NIW.marginals(ground, ground, "naive") #DEBUG

    message("Starting test runs:")
    
    for(alg in c("snappy1", "snappy2", "snappy3")) {
        message("Beginning ",alg,"...")
        try({ #wrap so that one algorithm being broken doesn't break the others (but we'll still see error messages from them)
        
            rNIW = get(paste("rNIW", alg, sep="."))
            
            tic = proc.time()
            samples = rNIW(n, kMu, kKappa, kPsi, kDF)
            toc = proc.time()
            runtime = (toc - tic)["elapsed"] #XXX this is an indulgence; runtimes should be done systematically with test.speed.main.R
            message("Runtime for ", alg,"(",n,",...): ", runtime,"s")
            
            # make density plots
            plot.NIW.marginals(ground, samples, alg)
            
            # make first moment convergence plots
            # TODO ....
            
            # make second moment convergence plots
            # TODO ....
            
            
            # test for equidistribution numerically 
            ks.test.NIW.marginals(ground, samples, alg)
            
            
        })
        
        message()
    }
    message("------------------------------------")
}


# our tests of interest:
# 1) sample values
#  -> marginal densities
#  -> moments
# 2) runtimes
# these purposes feel like they should share code
# but do not be overzealous: 
# runtimes 

main()
