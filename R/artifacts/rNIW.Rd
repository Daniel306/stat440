\name{rNIW}
\alias{rNIW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
Package has functions relating to the Normal-Inverse Wishart distribution.
}
\usage{
rNIW(n, Mu, Kappa, Psi, df)
%%rMNIW(n, Mu, Kappa, Psi, Df)
%%dNIW()
%%dMNIW()
%%lm.multivariable()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    the desired number of sample "points" (each "point" is a pair of a d-length vector and a dxd matrix)
  }
  \item{Mu}{
     a vector giving the mean of X samples
  }
  \item{Kappa}{
    a scalar which adjusts the variances of X
  }
  \item{Psi}{
  %% this one seems of
    a matrix, the mean of V samples
  }
  \item{df}{
    the "degrees of freedom" of the underlying Wishart distribution
  } 
  % Should probably also write down conditions on dimensions, maybe in details
}
\details{
 The dimensions match: d := length(Mu) == dim(Psi)[1] == dim(Psi)[2];
 
       Kappa is positive.

Psi is is positive-definite (=> square, symmetric, all positive eigenvalues, and chol(Psi) exists).

df need not be integer, but it must be larger than d-1.
}
\value{

%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{$X}{c(d, n) matrix of sample X values.}
  \item{$V}{c(d, d, n) matrix of sample V values.}
%% ...
}
\references{
http://en.wikipedia.org/wiki/Normal-inverse-Wishart_distribution
}
\author{
Daniel Galperin and Nick Guenther
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
rNIW is a general case of rMNIW %% add link to it
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

Psi <- matrix(c(2,1,1,2),2,2)
rNIW (n=10, Mu=c(1,2), Kappa = 1, Psi, df = 10) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
