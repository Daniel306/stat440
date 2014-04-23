setwd("~/GitHub/stat440/R/Rcpp")
require("Rcpp")
Rcpp::sourceCpp("test.cpp") #compile test.cpp on the fly and include it

a <- matrix(1:12,3,4)
t(a) %*% a
Mmult(t(a),a)
