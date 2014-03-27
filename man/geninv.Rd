\name{geninv}
\alias{geninv}
\title{Compute Moore-Penrose generalized inverse of a matrix}
\usage{
geninv(A)
}
\arguments{
  \item{A}{matrix}
}
\value{
The pseudoinverse of matrix A.
}
\description{
Compute Moore-Penrose generalized inverse of a matrix
}
\examples{
n <- 1000
p <- 500

x <- matrix(rnorm(n * (p - 1)), n, p)
x <- cbind(x, rowMeans(x))

## compute X'X
xpx <- crossprodcpp(x)

## compute generalized inverse of X'X
inv <- geninv(xpx)

## check if we have computed the M-P generalized inverse
all.equal(xpx, xpx \%*\% inv \%*\% xpx)
}

