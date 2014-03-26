Normal-Inverse-Wishart Generating Project
=========================================

Goal
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


Documentation
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