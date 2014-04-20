#File for testing rNIW.cpp

setwd("~/GitHub/stat440/R")

#source("NIW.R")
require("Rcpp")
Rcpp::sourceCpp("rNIW.cpp")

n = 10;
Mu= c(1,2);
Kappa= 1;
Psi=matrix(c(1,0,0,1), 2, 2);
df= 2;
rNIW.Rcpp2(n, Mu, Kappa, Psi, df)
#rNIW_Rcpp_2(2,2,4,4,5,6)

# Test for backSolveInverse
#Rcpp::sourceCpp("rNIW.cpp")
#A <- BartlettFactorCpp(3,10);
#A.inv <- backSolveInverse(A,3);
#A.inv %*% A
