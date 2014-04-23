setwd("~/GitHub/stat440/R/Rcpp")
require("Rcpp")
Rcpp::sourceCpp("test.cpp") #compile test.cpp on the fly and include it

a <- matrix(1:9,3,3)
t(a) %*% a
Mmult(t(a),a)
