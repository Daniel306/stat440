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
            "version1",
            "version2",
            "version3",
            "Rcpp",
            "extremelynaive.RcppEigen",
            "version.RcppEigen"
            )

krNIWAlgorithms.colors = c(
          "grey", "black",
          "red", "red", "red",
          "blue", "blue", "blue"
          )

krNIWAlgorithms.linestyle = c(
          "dashed", "solid",
          "solid", "dotted", "dashed",
          "solid", "dotted", "dashed"
          )

# sketchily merge the three vectors above into a data.frame to make it easier to deal with
# there's a bug
krNIWAlgorithms = data.frame(algorithm=krNIWAlgorithms, color=krNIWAlgorithms.colors, linestyle=krNIWAlgorithms.linestyle)
rm(krNIWAlgorithms.colors, krNIWAlgorithms.linestyle)



#kS = rWishart(69, kDF, kPsi) #a high dimensional matrix for checking syntax on