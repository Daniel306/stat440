


http://scythe.wustl.edu/ <-- " Programs written using Scythe are generally much faster than those written in commonly used interpreted languages, 
http://gallery.rcpp.org/articles/dmvnorm_arma/ 

[WP:Invertible Matrix](https://en.wikipedia.org/wiki/Invertible_matrix#Blockwise_inversion) has a surprising amount of insight.

[Tim Davis](http://www.cise.ufl.edu/research/sparse/), whose sparse matrix solvers underly everything that's fast in numerical computation: matlab, scipy, CUDA, and R.

Dahtah on [Why an inverse-Wishart prior may not be such a good idea](http://dahtah.wordpress.com/2012/03/07/why-an-inverse-wishart-prior-may-not-be-such-a-good-idea/).  
Datah's [follow up](https://dahtah.wordpress.com/2012/08/22/priors-of-convenience/).  
[and](http://andrewgelman.com/2012/08/29/more-on-scaled-inverse-wishart-and-prior-independence/) some [others](http://www.themattsimpson.com/2012/08/20/prior-distributions-for-covariance-matrices-the-scaled-inverse-wishart-prior/)  

[AND](http://www.statsblogs.com/2012/08/22/the-scaled-inverse-wishart-prior-distribution-for-a-covariance-matrix-in-a-hierarchical-model/)

http://mathoverflow.net/questions/160206/conditional-distribution-of-inverse-wishart

[andrew gelman's semi related thingy](http://andrewgelman.com/2006/09/01/modeling_the_gr/)

http://stackoverflow.com/questions/18349053/fastest-way-for-multiplying-a-matrix-to-a-vector:

* mmult (??)
* sweep (??)
* %*%
* t(t(mat)*v))

[IBM's BLAS docs](http://publib.boulder.ibm.com/infocenter/clresctr/vxrx/index.jsp?topic=%2Fcom.ibm.cluster.essl.v5r2.essl100.doc%2Fam5gr_sksubs.htm)



MCMC
------

* [in several languages](http://darrenjw.wordpress.com/2011/07/16/gibbs-sampler-in-various-languages-revisited/)
* [in cython](http://pyinsci.blogspot.ca/2010/12/efficcient-mcmc-in-python.html), merging elegant code and fast runtimes
