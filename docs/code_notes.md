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
