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
    #n = 33 #DEBUG

    m = 15500
    #m = 15 #DEBUG
    
    message()
    message("------------------------------------")
    message("Generating reference samples.")
    ground = rNIW.naive(m, kMu, kKappa, kPsi, kDF)
    #plot.NIW.marginals(ground, ground, alg="naive") #DEBUG
    #ks.test.NIW.marginals(ground, ground, "naive") #DEBUG

    message("Starting test runs:")
    
    for(alg in krNIWAlgorithms) { 
        if(alg == "extremelynaive") next; #hack: this one is REALLY SLOW, and only useful as a benchmark, so we drop it explicitly here.
        if(alg == "naive") next;          #hack: and this one is our reference samples, so its sort of silly to compare it.
        
        message("Beginning ",alg,"...")
        try({ #wrap so that one algorithm being broken doesn't break the others (but we'll still see error messages from them)
        
            rNIW = get(paste("rNIW", alg, sep="."))
            
            tic = proc.time()
            samples = rNIW(n, kMu, kKappa, kPsi, kDF)
            toc = proc.time()
            runtime = (toc - tic)["elapsed"] #XXX this is an indulgence; runtimes should be done systematically with test.speed.main.R
            message("Runtime for ", alg,"(",n,",...): ", runtime,"s")

            ## (marginal) DISTRIBUTIONS
            message("Density plots")
            # 1) computationally
            plot.densities(ground$X, samples$X, "NIW X")
            plot.densities(ground$V, samples$V, "NIW V")
            # 2) analytically
            #analytic_densities = NIW.densities(kMu, kKappa, kPsi, kDF)
            # ^ does this step return like.. a.. tuple of matrices of functions?
            # does R even support such a beast?
            # ..I might have to move the marginals_do to this function
            #plot.densities(analytic_densities$X, samples$X) #TODO: work out just how this is...
            #plot.densities(analytic_densities$V, samples$V) #I only know the diagonal not know what this is, so skip it
            
            
            ## MOMENTS
            # make first moment convergence plots
            message("Moments")
            # first moments (matrix and non-matrix)
            # 1) computationally
            #message("m1comp") #DEBUG
            #plot.mean.convergence(ground$X, samples$X, "NIW X")
            #plot.mean.convergence(ground$V, samples$V, "NIW V")
            # 2) analytically
            #message("m1analytic") #DEBUG
            analytic = NIW.mean(kMu, kKappa, kPsi, kDF)
            plot.convergence(analytic$X, samples$X, mean, "NIW X", sub="(analytically)")
            plot.convergence(analytic$V, samples$V, mean, "NIW V", sub="(analytically)")
                        
            # first variances
            # 1) computationally
            #message("v1comp") #DEBUG
            plot.var.convergence(ground$X, samples$X, "NIW X")
            plot.var.convergence(ground$V, samples$V, "NIW V")
            # 2) analytically
            #message("v1analytic") #DEBUG
            analytic = NIW.var(kMu, kKappa, kPsi, kDF)
            plot.var.convergence(analytic$X, samples$X, "NIW X", sub="(analytically)")
            plot.var.convergence(analytic$V, samples$V, "NIW V", sub="(analytically)")
            #  ^ this specialcase is going to bite!
            
            # second matrix moments: the means of the outer products
            # 1) computationally
            #message("mm1comp") #DEBUG
            plot.mean.convergence(fancy_matrix_square(ground$X), fancy_matrix_square(samples$X), "XX^T")
            plot.var.convergence(fancy_matrix_square(ground$V), fancy_matrix_square(samples$V), "VV^T")
            # 2) we do not do this analytically
            # that's all she wrote, folks!
            
            # variances of the second matrix moments are unknown and uninteresting
            #  --the 2nd moments themselves are already close to a variance.
            # None of this handles element-element covariance, which is a whole other headtrip.
            
            ## KS TESTS
            # this tests for equidistribution numerically
            message("KS test results for ", alg) #this is sort of awkward, that there's message()s in this call
            # TODO: rewrite using marginals_do
            NIW.ks.test(ground, samples)
                        
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
