\name{dMNIW}
\alias{dMNIW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DENSITY OF MNIW
}
\description{
calculates density of MNIW based on X, V, and the parameters associated with its distribution.
}
\usage{
dNIW(X, V, Mu, Kappa, Psi, df, log=f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X,V}{ points of vectors to test}
  \item{Mu,Kappa,Psi,df}{parameters of NIW}
  \item{log}{whether to produce log-probabilities or not}
}

}
\details{
X is (d,q) matrix and V is a (d,d) matrix
or X is a (d,q,n) matrix and V is a (d,q,n) matrix
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
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
X <- Matrix(c(1,2,3,4), 2,2)
V <- matrix(c(3,1,1,3), 2,2)
Psi <- diag(2)
Mu <- matrix(0,2,2)
Kappa <- diag(2)
df <- 10
dMNIW(X,V, Mu = Mu, Kappa = Kappa, Psi = Psi, df = df)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line