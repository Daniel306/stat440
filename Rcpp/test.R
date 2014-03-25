require("Rcpp")
Rcpp::sourceCpp("test.cpp") #compile test.cpp on the fly and include it

foo()
