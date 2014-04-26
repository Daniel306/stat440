Normal-Inverse-Wishart Generating Project
=========================================

Goals
----

The [NIW](https://en.wikipedia.org/wiki/Normal-inverse-Wishart_distribution) distribution
is a typical conjugate prior for the multivariate normal distribution,
and using in R is currently very, very, very slow.

Here's a real working statistician [running up against it in R](https://dahtah.wordpress.com/2012/03/07/why-an-inverse-wishart-prior-may-not-be-such-a-good-idea/):
```{R}
 #Generate n samples from the prior
 rprior <- function(n=1,r=3,M=diag(rep(1,2)))
 {
 Minv <- solve(M)
 rlply(n,chol2inv(chol(rwish(r,Minv))))
 }
 
#Wishart samples
 rwish <- function(r,R)
 {
 X <- rmvnorm(r,sig=R)
 t(X)%*%X
 }
```

The above algorithm is terribly inefficient: it literally constructs Wishart matrices from their definition and inverting them n times.
That snippet does not even include the part where she samples the corresponding Normals,
which requires refactoring the variance matrix (the `chol` above).

This package does better by exploiting a) algebraic shortcuts b) [Rcpp](http://rcpp.org).

We are aware of two other Inverse-Wishart implementations for R,
`riwish` of [MCMCpack](http://mcmcpack.wustl.edu/) and 
[`rinvwishart`](https://github.com/Statisticat/LaplacesDemon/blob/master/R/distributions.R)
 of [LaplacesDemon](http://www.bayesian-inference.com/software),
 and we are also aware of the standard Wishart implementation
 [`rWishart`](http://svn.r-project.org/R/trunk/src/library/stats/src/rWishart.c()) in base R.
 This package is different because it is fast enough to do heavy bayesian
 statistics in the environment statisticians are already using.

To make its use clear, this package also includes an application of the sampler to
[multivariable regression](https://en.wikipedia.org/wiki/Bayesian_multivariate_linear_regression).

<!--  [the](http://www.ats.ucla.edu/stat/stata/dae/mvreg.htm) [alternatives](http://cameron.econ.ucdavis.edu/excel/ex61multipleregression.html). -->

This is **alpha** code.

Code Guide
----------------

### Usage

This is alpha code because it is not packaged as an R package yet.

To use, make sure you have the [dependencies](#dependencies) installed and run the tests:
```
$ Rscript test.correctness.main.R
$ Rscript test.speed.main.R
```

Then copy the source to your working directory, and run
```{r}
source("dNIW.R")
source("rNIW.R")
source("MultivariableRegression.R")
```

This will give you access to these APIs:

Densities (pdfs):

* `dNIW(X, V, Mu, Kappa, Psi, df, log=F)`
* `dMNIW(X, V, Mu, Kappa, Psi, df, log=F)`

Samplers

* `rNIW(n, Mu, Kappa, Psi, df)`
* `rMNIW(n, Mu, Kappa, Psi, df)`
* `rmultivariableregression(points, B, V)`

Multivariable Regression

* `lm.multivariable(m, X, Y, Lambda, Omega, Psi, df)`

The MNIW functions are a further generalization of this generalized distribution
 operating with [Matrix-Normals](https://en.wikipedia.org/wiki/Matrix_normal_distribution).
 
The "regression" operates in a Bayesian way: it returns `m` samples of the estimated parameters `B` and `V`.
 
See the Rdocs (e.g. `help(rNIW)`) for all the details.

### Files

* `rNIW.R` contains the various sampler implementations we experimented with
* `dNIW.R` contains analytic density functions for the distribution
* `MultivariableRegression.R` contains the "multivariable regression"
* `test.*.R` is code to vet the implementations above.
    * `test.speed.R` is a harness to record runtimes of the different implementations, plot them, and summarize the results in a table.
    * `test.correctness.R`  vets that this code is actually sampling the proper distribution by reducing the problem to look at the marginals. It provides
        * _analytic_ and _computational_ (kernel density estimate) density functions.
        * _analytic_ and _computational_ (sample mean) moments.
        * _computational_ (KS-test) distribution comparison.
    * `test.MultivariableRegression.R` has some basic examples and smoketests of `lm.multivariable()`
    * `test.*.main.R` are scripts which can be run directly to exercise all of the above:
        * `test.speed.main.R` produces a plot and a table of runtimes
        * `test.correctness.main.R` produces a very large number of histogram-density plots and converging moment plots, and text streams of results from the KS-tests (by default at an alpha=0.05)
        * `test.MultivariableRegression.main.R` runs `lm.multivariable()` and prints the resulting point estimate on
            * a) artificial data with known parameters
            * b) the Boston housing dataset in MASS

### Dependencies

* [R](http://r-project.org)
* [Rcpp](http://cran.r-project.org/web/packages/Rcpp/index.html)
* [RcppEigen](http://cran.r-project.org/web/packages/RcppEigen/index.html)


Future Work
-----------

1. Parallelize the algorithm, and investigate on GPUs
1. Allow the optional cholesky parameterization, and use that to reduce the number of algorithms
1. Investgate if this method is faster than the naive implementation with a good `plyr` call.
1. Actually create a package.
1. Make lm.multivariable return a more useful structure.