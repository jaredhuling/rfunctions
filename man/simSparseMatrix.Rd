\name{simSparseMatrix}
\alias{simSparseMatrix}
\title{Simulate a matrix with very few nonzero entries}
\usage{
simSparseMatrix(sparsity, dim, boolean = T, slam = F)
}
\arguments{
  \item{sparsity}{numeric value between 0 and 1 indicating
  level of sparseness}

  \item{dim}{numeric vector of length 2 indicating matrix
  dimensions}

  \item{boolean}{TRUE/FALSE. if true, matrix will only take
  values 0 and 1}

  \item{slam}{if TRUE, will return matrix in
  simple_triplet_matrix form. If FALSE, will return
  SparseMatrix object}
}
\value{
SparseMatrix or simple_triplet_matrix
}
\description{
Simulate a matrix with very few nonzero entries
}
\examples{
n.obs <- 1e5
n.vars <- 1e3

#simulate a very sparse matrix (this matrix has many zeros and few ones)
x.sparse <- simSparseMatrix(sparsity = 0.99, dim = c(n.obs, n.vars), boolean = T)
}

