\name{dNIW}
\alias{dNIW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DENSITY OF NORMAL-INV-WISHART
}
\description{
calculates density of NIW based on X, V, and the parameters associated with its distribution.
}
\usage{
dNIW(X, V, Mu, Kappa, Psi, df, log=f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X,V}{ points of vectors to test}
  \item{Mu,Kappa,Psi,df}{paramaters of NIW}
  \item{log}{whether to produce log-probabilities or not}
}

}
\details{
X is d-length vector or a (d,1) matrix and V is a (d,d) matrix
or X is a (d,n) matrix and V is a (d,d,n) matrix
}
\value{
Returns a vector of probabilities p or log(p) if log==TRUE of length n
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Daniel Galperin and Nick Guenther
}
\note{
Related to rNIW and dMNIW
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
X <- c(1,2)
V <- matrix(c(3,1,1,3), 2,2)
Psi <- diag(2)
df <- 10
dNIW(X,V, Mu = c(1,2), Kappa = 1, Psi = Psi, df = df)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
