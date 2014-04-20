#File for testing rNIW.cpp

setwd("~/GitHub/stat440/R")

source("NIW.R")
source ("test.constants.R")
require("Rcpp")
Rcpp::sourceCpp("rNIW.cpp")
n=100000;
sample1 <- rNIW.Rcpp2(n, kMu, kKappa, kPsi, kDF)
#rNIW_Rcpp_2(n, 2, Mu, Kappa, Psi, df)

# Test for backSolveInverse
#Rcpp::sourceCpp("rNIW.cpp")
# A<- BartlettFactorCpp(2,10);
#A.inv <- backSolveInverse(A,3);
#A.inv %*% A
sample2 <-rNIW.snappy2(n, kMu, kKappa, kPsi, kDF)

rowMeans(sample1$X)
rowMeans(sample2$X)
