\name{crossprodcpp}
\alias{crossprodcpp}
\title{Compute X'X or X'WX where W is diagonal weight matrix}
\usage{
crossprodcpp(x, w = NULL)
}
\arguments{
  \item{x}{matrix input}

  \item{w}{numeric vector. If specified, X'diag(w)X will be
  computed}
}
\value{
A symmetric matrix X'X
}
\description{
Compute X'X or X'WX where W is diagonal weight matrix
}
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

