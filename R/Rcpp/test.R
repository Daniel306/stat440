setwd("~/GitHub/stat440/R/Rcpp")
require("Rcpp")
Rcpp::sourceCpp("test.cpp") #compile test.cpp on the fly and include it

a <- matrix(1:6,3,2)

foo(a)
