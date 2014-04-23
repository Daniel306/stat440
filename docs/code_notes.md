Scrap notes on codey things
=============================

Double-Sourcing Rcpp
-----------------------

R is clever enough about not reloading things unnecessarily that this code is speedy:

```
f = function() { 
  require("Rcpp")
  sourceCpp("rNIW.cpp")
 }
 
f()
# Loading required package: Rcpp
f() # <-- does not recompile rNIW.cpp
f() # <-- this neither
```


Using pointers in Rcpp
---------------------------



Inverting
---------------

The honest to god fastest way to invert a triangular matrix is (apparently, eg [see here](http://gallery.rcpp.org/articles/dmvnorm_arma/)) 
```{r}
backsolve(M, diag(length(diag(tri))))
```

and of course the fastest usual way is
```{r}
solve(M)
```

`diag` is quirky: given a vector, it makes a matrix with that as its diagonal; given a matrix it extracts the diagonal and returns it as a vector (this mirrors matlab behaviour). But given a *scalar*, it creates the identity matrix of that size.

