\name{rmultivariableregression}
\alias{rmultivariableregression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generates X and Y
}
\description{
 generates from is Y = XB + E where (E_i)^T ~iid MultiNormal(0, V)
}
\usage{
rmultivariableregression(points, B,V)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{points}{Either an (n,p) matrix X, or or a scalar n}
  \item{B}{B values of model. Should be c(p,q) }
  \item{B}{Y values of model. Should be c(q,n) }
}
\details{
n is either the inter input for points or the
}
\value{

  \item{X}{c(p,n) matrix }
  \item{Y}{c(q,n) matrix }
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
lmMultivariate() which generates B,V from X,Y and priors
}
\examples{
rmultivariableregression(5, diag(1), diag(1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
