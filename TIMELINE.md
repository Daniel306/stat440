Normal-Inverse-Wishart Generating Project: TODO and Timeline
============================================================

D-7 (Wed)
---------

[daniel] [nick] sync


D-8 (Thur)
----------

[nick]

- [x] freeze rNIW API
    - [x] figure out if we should be c(n,d,d) or c(d,d,n) 
    - [x] make the typecheck wrapping thing..better
    - [x] update implementations to match frozen API
- [x] testing harness for plotting marginals (req: generating ground truth)
- [x] testing harness for ks tests
- [x] test each current impl with the tests
- [x] REPAIR each current impl we know all the ones that assume chol(inv(P)) = inv(chol(P)) are wrong.
    - [ ] save the wrong ones as examples to paste into the report 

[daniel]

- [x] read up on the examples assigned to us


D-7 (Friday)
------------

[daniel] [nick]

- [-] meet 12:30s

[nick]

- [x] runtime testing harness
    - [x] measure over n
    - [ ] measure over d (requires randomization to be in place) 
    

[daniel]

- [x] testing harness plotting moments (specifically, convergence of moments)

D-6 (sat)
----------

[nick]  

- [x] tighten the moment tests
    - [x] simplify the code to use a single function mean.convergence(M, dim) which takes cumulative means over dim 'dim' (probably via `apply`)
    - [x] use crossprod() instead of %*%; use apply(cumsum) instead of cumsum in a loop, etc
    - [x] use lines() instead of plot bc plotting 1e5 points is absurd and is killing our machines to render
        - instead, chose to write intseq() to reduce the number of plotted points
    - [x] plot the expected founds from ground
        - [x] special case ground so that a SINGLE value is taken to be the expected value 
- [x] repair what `mathematics.tex` thinks about chol()

[daniel]

- [x] have NIW Rcpp impl written

D-5 (Sunday)
----

[nick] 

- [x] read (and struggle with) Lysy's answer 1

[daniel]

- [ ] rename implementations from 'snappy' to something more professional

D-4 (Monday)
---

[nick]

- [ ] randomize the parameters used for testing (requires test.correctness to be able to handle arbitrary ground truths..)
    - [ ] random positive definite matrix' part, which would be good to write as a subroutine...)
- [ xish] dNIW() 
    - [x] freeze dNIW API
    - [x] implementation
    - [ ] come up with a testing plan for dNIW
        - easiest: take *conditionals* (which are *proportional* to the marginals?? <-- no, this is NOT true)
        - harder: ..take marginals out of it? this requires using integrate() over multiple dimensions to get the correct scaling factor
        - [x] convince myself that the conditional distributions are different from the marginals, and there's no computationally cheap way to convert one to the other
- [ ] read Lysy's answer 2
- [x] repair my R CRAN library after I updated and went 3.0->3.1
- [x] write marginalize()
    - [ ] write NIW.X.marginal
    - [ ] write NIW.V.marginal ----- ooo this one is tricky because V has to be symmetric
    - [ ] use this to do dNIW testing
    ^ PROBLEM: marginalize is absurdly slow; it's useless for more than 3 covariates

D-3 (Tuesday)
----

[nick]

- [x] MNIW.typecheck()
    - [ ] factor NIW.typecheck so that it shares code with that
- [x] split NIW.R to dNIW.R, dMNIW.R, rNIW.R and rMNIW.R
- [x] clean up this TODO list to git
- [x] write naive version in Rcpp (for timing comparison)
- [ ] in the runtime tests, map the runtime matrix to have a ratio column so that we can say "snappy2 is 4.3 times faster than naive"
- [ ] factor the algorithm list to `test.constants.R`
- [ ] figure out if using pointers/references in Rcpp is faster; specifically, does passing a NumericVector cause a COPY of that vector even in C?

[daniel]

- [x] (9am-11:30) ACTSC exam
- [ ] rmNIW
- [ ] dmNIW (matrix Normal inverse Wishart)
- 
- [ ] i.i.d. Matrix-Normal Sampler
- [ ] Multivariable Regression algorithm
- [ ] Test multivar regression against the Matrix-Normal sampler

D-2 (Wednesday)
---

[nick]
- [ ] Linear Algebra parlour tricks:
    - [ ] factor the common matrix terms to before/after the loop (call this `snappy4`)
    - [ ] does this save any time??
- [ ] 13:00 - 14:00 blocked off
- [ ] Dig up analytic formula for the particular marginals of NIW and (note: there's 2^(# parameters) different marginals; pick wisely which to look at)
    - [ ] implement as functions
    - [ ] wrap such functions into thingies and use them as 'ground' for the marginal plots
    - t-distributions (and the matrix-normal has multi-t)
        - analytically
        - [x] against ground samples

[daniel]

- [ ] double-check (or possibly, first-check) the marginals nick dug up

[daniel] [nick]

- [ ] Mathematics writeup vetted and tightened; relevant parts clipped into `report.tex`

D-1 (Thursday)
----

[daniel]

- [ ] Request sample real world dataset.
- [ ] Apply Multivariable Regression to sample real world dataset.

[nick]

- [ ] make R package(s) (i.e. call require() instead of Rcpp::source())
    - [ ] NIW
        - [ ] dNIW()
        - [ ] rNIW()
        - [ ] dMNIW()
        - [ ] rMNIW()
        - +various utility functions ((maybe these should get their own package?))
    - [ ] MultivariableRegression
        - depends: NIW
        - [ ] ???
    - ([ ] HierarchicalMatrixNormalGibbsSampler  _there is no way we're doing this by Friday_)
    - [ ] prototype + testing
    - [ ] write R help files for each function
- [ ] extract the prototype code (i.e., apart from the final package), package it up reusably

D (Friday, April the 25th, 2014)
---------------------------------

[daniel] [nick]

- [ ] completed `report.tex`:
    - [ ] code (or at least the relevant parts) prettyprinted and included in the appendix (use knitr for this??)
    - [ ] zipfile with the remainder of the code 
- [ ] upload report to Learn or email'd to Lysy or whatever

D+k
---------------------

- [ ] move to rbenchmark::benchmark instead of NIW.runtimes
- [ ] write vignette()s
- [ ] look into the curve() function, instead of lines()
- [ ] derive how the NIW is a conjugate prior for a hierarchical normal
    - [ ] record in `mathematics.tex`
- [ ] Hierarchical normal regression via Gibbs sampling (this is a separate piece of code; it may share some subroutines, but its core is an entirely different beast)
- [ ] write versions of every function which are cholesky-parameterized (credit to LaplacesDemon for the simple and elegant idea);
    - [ ] use these versions to collapse large portions of the code
- [ ] move to using RUnit
- [ ] test for and quantify numerical instability
- [ ] runtimes
    - [ ] measure over d (requires randomization to be in place) 

Open Questions
================
- how do we test dNIW? `marginalize()` is too slow to be used on more than 3 variates, which outruns the dimensionality of the smallest interesting NIW distribution
- [x] the existence of the kroenecker product breaks the rNIW and dNIW APIs
     - **Answer**: yes, it does; change the API: rMNIW() for the Matrix-Normal|Inverse-Wishart distribution 
- [x] what is the proper IW PDF? is wikipedia wrong?
    - **Answer**: probably not; recall that changes of variables necessarily involve a scaling factor--the Jacobian--which pulls in a determinant; the particular algebra is beyond our ken at this point, though.
- [ ] the iid multivariate normals all share a single variance matrix V; please re-explain what this means; shouldn't the variances be scaled? What is the interpretation of this model?
- [ ] Why is snappy2 faster than snappy3? Shouldn't using backsolve()s be faster than %*%?
- [ ] Is psi the precision matrix or the variance matrix? Why do you invert a matrix to get a covariance-esque matrix? What is going on?
- [ ] "lint" the code:
    - [ ] naming conventions are followed
        - [ ] distribution parameters are capitalized, except for df
        - [ ] density functions begin with 'd'
        - [ ] samplers begin with 'r'
        - [ ] 
    - [ ] indentation in R is 2 spaces
    - [ ] indentation in C is 4 spaces
    - [ ] brace style (?)
    - [ ]  ????
    - [ ] find all the TODOs and deal with them

