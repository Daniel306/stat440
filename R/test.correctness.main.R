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

            ## DISTRIBUTIONS
            message("Marginal densities")
            # 1) computationally
            #plot.densities(ground$X, samples$X, "NIW X", sub=paste(alg))
            #plot.densities(ground$V, samples$V, "NIW V", sub=paste(alg))
            # 2) analytically
            #analytic_densities = NIW.densities(kMu, kKappa, kPsi, kDF)
            # ^ does this step return like.. a.. tuple of matrices of functions?
            # does R even support such a beast?
            # ..I might have to move the marginals_do to this function
            analytic = NIW.densities(kMu, kKappa, kPsi, kDF, n)
            plot.densities(analytic$X, samples$X, "NIW X", sub=paste(alg, "(analytically)"))
            plot.densities(analytic$V, samples$V, "NIW V", sub=paste(alg, "(analytically)"))
            
            ## MOMENTS
            # make first moment convergence plots
            message("Moments")
            # first moments (matrix and non-matrix)
            # 1) computationally
            #message("m1comp") #DEBUG
            #plot.means(ground$X, samples$X, "NIW X", sub=paste(alg))
            #plot.means(ground$V, samples$V, "NIW V", sub=paste(alg))
            # 2) analytically
            #message("m1analytic") #DEBUG
            analytic = NIW.mean(kMu, kKappa, kPsi, kDF)
            plot.means(analytic$X, samples$X, "NIW X", sub=paste(alg, "(analytically)"))
            plot.means(analytic$V, samples$V, "NIW V", sub=paste(alg, "(analytically)"))
              
            # first variances
            # 1) computationally
            #message("v1comp") #DEBUG
            #plot.vars(ground$X, samples$X, "NIW X", sub=paste(alg))
            #plot.vars(ground$V, samples$V, "NIW V", sub=paste(alg))
            # 2) analytically
            #message("v1analytic") #DEBU
            analytic = NIW.var(kMu, kKappa, kPsi, kDF)
            plot.vars(analytic$X, samples$X, "NIW X", sub=paste(alg, "(analytically)")) 
            plot.vars(analytic$V, samples$V, "NIW V", sub=paste(alg, "(analytically)"))
            #  ^ this specialcase is going to bite!
            
            # second matrix moments: the means of the outer products
            # 1) computationally
            #message("mm1comp") #DEBUG
            #plot.means(fancy_matrix_square(ground$X), fancy_matrix_square(samples$X), "XX^T", sub=paste(alg))
            #plot.means(fancy_matrix_square(ground$V), fancy_matrix_square(samples$V), "VV^T", sub=paste(alg))
            # 2) we do not do this analytically
            # [[that's all she wrote, folks!]]
            
            # variances of the second matrix moments are unknown and uninteresting
            #  --the 2nd moments themselves are already close to a variance.
            # 1) computationally
            # [[ see above ]]
            # 2) analytically
            # [[ see above ]]
            
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
