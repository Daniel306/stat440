# test.correctness.main.R
#
#  Until we learn the proper Rish way of writing tests.
#  This script will serve.
#  in Python, the contents of this file would be under the
#   "if __name__ == '__main__'" 
#  block.

source("test.correctness.R")
source("test.constants.R")

main <- function() { # avoid polluting the namespace
    
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
    
    for(alg in c("snappy1", "snappy2", "snappy3", "Rcpp2")) { 
        message("Beginning ",alg,"...")
        try({ #wrap so that one algorithm being broken doesn't break the others (but we'll still see error messages from them)
        
            rNIW = get(paste("rNIW", alg, sep="."))
            
            tic = proc.time()
            samples = rNIW(n, kMu, kKappa, kPsi, kDF)
            toc = proc.time()
            runtime = (toc - tic)["elapsed"] #XXX this is an indulgence; runtimes should be done systematically with test.speed.main.R
            message("Runtime for ", alg,"(",n,",...): ", runtime,"s")
            
            # make density plots
            plot.NIW.densities(ground, samples, alg)
            
            # make first moment convergence plots
            plot.NIW.moment.first(ground, samples)
                        
            # make second moment convergence plots
            plot.NIW.moment.second(ground, samples)
            
            # test for equidistribution numerically 
            NIW.ks.test(ground, samples, alg)
            
            
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
