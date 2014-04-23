#File for testing rNIW.cpp

setwd("~/GitHub/stat440/R")

source("NIW.R")
source ("test.constants.R")
require("Rcpp")
Rcpp::sourceCpp("rNIW.cpp")

# test means:
#n=2000000;
sample1 <- rNIW.Rcpp2(n, kMu, kKappa, kPsi, kDF);
#rowMeans(sample1$X)
#sample2 <- rNIW.snappy2(n, kMu, kKappa, kPsi, kDF);
#rowMeans(sample2$X)
# Appears correct:
# > rowMeans(sample1$X)
#[1]  0.2461772 -1.5002524  0.2891439  8.9968747
#> sample2 <- rNIW.snappy2(n, kMu, kKappa, kPsi, kDF);
#> rowMeans(sample2$X)
#[1]  0.2514245 -1.4995739  0.3229725  8.9993225


# tested the multiplication of inverses (by extracting values and doing calculation in R too)
rNIW.Rcpp2(1, kMu, kKappa, kPsi, kDF);
rNIW.snappy2(1, kMu, kKappa, kPsi, kDF);
# Test for bartlett decomposition
#Rcpp::sourceCpp("rNIW.cpp")
#sample1 <- matrix(NA, 4*4*n);
#dim(sample1) <- c(4,4,n);
#sample2 <- matrix(NA, 4*4*n);
#dim(sample2) <- c(4,4,n);

#for(i in 1:n){
 # sample1[,,i] = BartlettFactorCpp(4,kDF);
  #sample2[,,i] = BartlettFactor(4,kDF);  
#}
#rowSums(sample1, dims = 2)/n -rowSums(sample2, dims = 2)/n
