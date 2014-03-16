\name{gklBidiag}
\alias{gklBidiag}
\title{Compute Largest Singular Value of Matrix x}
\usage{
gklBidiag(x, v, maxit, upper.bound.prob = NULL)
}
\arguments{
  \item{x}{matrix input}

  \item{v}{numeric vector. Initialize GKL bidiagonalization
  with random vector on unit sphere}

  \item{maxit}{integer. Maximum number of Lanczos
  iterations}

  \item{upper.bound.prob}{upper bound probability for the
  largest singular value}
}
\value{
list
}
\description{
Compute Largest Singular Value of Matrix x
}
\examples{
n.obs <- 1e5
n.vars <- 150

x <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)

## compute largest singular value of x
lanczos <- gklBidiag(x, runif(ncol(x)), 10L, 0.99)
str(lanczos)
}

