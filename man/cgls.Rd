% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cgls.R
\name{cgls}
\alias{cgls}
\title{Solve the linear system (A'A + lambda * I)x = Ab using cgls}
\usage{
cgls(A, b, lambda = 0, maxit = 50L, tol = 1e-05, init = if
  (is.null(dim(b))) {     numeric(length(b)) } else {     array(0, dim = dim(b))
  })
}
\arguments{
\item{A}{matrix or CsparseMatrix input}

\item{b}{numeric vector or numeric matrix (corresponding to multiple right-hand sides).}

\item{lambda}{numeric positive scalar for ridge penalty}

\item{maxit}{integer. Maximum number of iterations}

\item{tol}{tolerance level}

\item{init}{initialization of estimate}
}
\value{
list
}
\description{
Solve the linear system (A'A + lambda * I)x = Ab using cgls
}
\examples{
n <- 10000
p <- 100

set.seed(100)
x <- simSparseMatrix(0.8, boolean = FALSE, dim = c(n, p)) 
b <- as.vector(x \%*\% rnorm(p) + rnorm(n))

A <- crossprod(x)
xtb <- as.vector(crossprod(x, b))

system.time( alpha.true <- solve(A, b) )
system.time( alpha.cgls <- cgls(x, b) )

max(abs(alpha.cgls$x - alpha.true))
}

