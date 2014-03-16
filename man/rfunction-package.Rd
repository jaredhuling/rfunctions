\name{rfunctions-package}
\alias{rfunctions-package}
\alias{rfunctions}
\docType{package}
\title{
A set of R functions that I find useful
}
\description{
R functions that are useful or are fast. Most of the implementation is in C++, utilizing RcppEigen
}
\details{
\tabular{ll}{
Package: \tab rfunctions\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-03-15\cr
License: \tab  GPL (>= 2)\cr
Depends: \tab Rcpp, RcppEigen\cr
}
oemfitFast is the main OEM function. It takes x and y as arguments for data and produces the regularization path over a grid of values for lambda. It accepts x as either a 'matrix' or 'sparseMatrix'
}
\author{
Jared Huling

Maintainer: Jared Huling <jaredhuling [ at ] gmail [ dot ] com>

}


\keyword{ package }

\examples{
n.obs <- 1e5
n.vars <- 150

x <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)

## compute X'X
xpx <- crossprodcpp(x)

weights <- runif(n.obs)

## compute X'WX
xpwx <- crossprodcpp(x, weights)


}
