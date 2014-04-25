Normal-Inverse-Wishart Generating Project
=========================================

Goals
----

The [NIW](https://en.wikipedia.org/wiki/Normal-inverse-Wishart_distribution) distribution is a typical, and is currently very, very, very slow. Here's a real working statistician [running up against it, painfully, in R](https://dahtah.wordpress.com/2012/03/07/why-an-inverse-wishart-prior-may-not-be-such-a-good-idea/):
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

We are aware of two other Normal|Inverse-Wishart implementations for R,
and the `[rWishart](http://svn.r-project.org/R/trunk/src/library/stats/src/rWishart.c())`
function in base R. This package is different because it is fast enough to
do heavy bayesian statistics in the environment statisticians are already using.

To be useful, this package also includes an application the sampler to
[multivariable regression](https://en.wikipedia.org/wiki/Bayesian_multivariate_linear_regression).

<!--  [the](http://www.ats.ucla.edu/stat/stata/dae/mvreg.htm) [alternatives](http://cameron.econ.ucdavis.edu/excel/ex61multipleregression.html). -->


Code Guide
----------------

### Files

* `rNIW.R` contains the various sampler implementations we tried
* `dNIW.R` contains analytic density functions for this distribution
* `MultivariableRegression.R` contains the a "multivariable regression"
* `test.*.R` is code to vet the implementations above.
    * `test.speed.R` is a harness to record runtimes of the different implementations, plot them, and summarize the results in a table.
    * `test.correctness.R`  vets that this code is actually sampling the proper distribution by reducing the problem to look at the marginals. It provides
        * _analytic_ and _computational_ (kernel density estimate) density functions
        * _analytic_ and _computational_ (sample mean) moments
        * _computational_ (KS-test) distribution comparison
                * `test.MultivariableRegression.R` 
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


How To Keep Documentation
--------------

`.md` files (like this one) are useful for documenting our work in a simple, **portable** format.
Here are some options for working with markdown:

1. Learn the syntax and using a plain text editor.
2. Directly on GitHub (_Web_)
3. [Lightpaper](http://clockworkengine.com/lightpaper-mac/) (_OS X / Android_)
3. [Texts](http://www.texts.io/) (_OS X_, _Windows_, with Math support)
3. [Mou](http://mouapp.com/) (_OS X_)s
4. [MarkdownPad](http://www.markdownpad.com/) (_Windows_)
5. [WriteMonkey](http://writemonkey.com/) (_Windows_)
6. [StackEdit](http://stackedit.io/) (_Web_; also supports inlining LaTeX trivially with `$ $`)
7. [reText](http://sourceforge.net/p/retext/home/ReText/) (_Linux_, _OS X_)


$$ x=\frac{-b \pm \sqrt {b^2-4ac}}{2a} $$

Keeping notes as we go will be vital. Keeping track of even the small things can save us. Because we're using git, remembering all of this, searching it, erasing it and then reviving it are all cheap and fearless.



Future Work
-----------

1. Extensions to GPUs
1. Allow the optional cholesky parameterization, and use that to reduce the number of algorithms
1. Investgate if this method is faster than the naive implementation with a good `plyr` call.
1. Actually create a package.
