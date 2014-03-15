\name{crossprodChunk}
\alias{crossprodChunk}
\title{Compute X'X in chunks}
\usage{
crossprodChunk(X, row.chunk = NULL, sparse = FALSE)
}
\arguments{
  \item{X}{matrix input}

  \item{row.chunk}{integer. Number of rows per chunk.
  Specify NULL for no chunking}

  \item{sparse}{Boolean value to specify if matrix X is
  sparse}
}
\value{
A symmetric matrix X'X
}
\description{
Compute X'X in chunks
}
\examples{
crossprodChunk(matrix(rnorm(n.obs * n.vars), 1e4, 1e2), row.chunk = 500, FALSE)
}

