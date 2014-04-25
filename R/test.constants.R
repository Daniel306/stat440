#test.constants.R
# particular Normal-Inverse-Wishart distribution parameters
# for testing

# TODO: let these be user set and/or randomize them
kMu = c(0.25, -1.5, 0.33, 9)
kKappa = 1
kPsi  = cbind(c(2.84, 0.43, 0.16, .5), c(0.43, 1.52, -0.24, -.73), c(.16, -.24, 4.49, 88), c(1,-2,3,4))
kPsi  = t(kPsi)%*%kPsi #after I added digits (to go from p=3 to p=4) at random,
                       #coerce the matrix to positive definite again, in the most crude way I know
kDF = 7.32 # must be at least as large as the dimension of kV

krNIWAlgorithms = c(
            "extremelynaive",
            "naive",
            "snappy1",
            "snappy2",
            "snappy3",
            "Rcpp2",
            "extremelynaive.RcppEigen",
            "snappy.RcppEigen"
            )

#kS = rWishart(69, kDF, kPsi) #a high dimensional matrix for checking syntax on