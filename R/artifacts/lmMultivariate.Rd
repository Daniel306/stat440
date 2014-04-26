\name{lm.multivariable}
\alias{lm.multivariable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian samples of Multivariate regression
}
\description{
generates Bayseian samples for B and V in a multivariate regression using rMNIW.
}
\usage{
lm.multivariable(m, X, Y, Lambda=NULL, 
                 Omega=NULL, Psi=NULL, df=71)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{ number of samples to draw}
  \item{X}{X values of model. Should be c(p,n) }
  \item{Y}{Y values of model. Should be c(q,n) }
  \item{Lambda}{prior mean of B, c(p,q)}
  \item{Omega}{Scale factor of variance, c(p,p)}
  \item{Psi}{prior mean of V, c(q,q)}
  \item{df}{degrees of freedom}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  \item{B}{c(p,q,n) matrix containing samples for B}
  \item{V}{c(q,q,n) matrix containing samples for V}
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
rMNIW which this calls
rmultivariableregression() which generates X,Y from B,V
}
\examples{
require(MASS)
D <- Boston
Y <- D[c("crim", "zn", "indus", "tax", "black", "lstat")]
X <- D[c("chas", "nox", "rm", "age", 
        "dis", "rad", "ptratio", "medv")]

# One example using flat priors and df of 73

Result <- lm.multivariable(100000,as.matrix(X),as.matrix(Y))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
