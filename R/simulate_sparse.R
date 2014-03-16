
#' Compute Largest Singular Value of Matrix x
#'
#' @param sparsity numeric value between 0 and 1 indicating level of sparseness
#' @param dim numeric vector of length 2 indicating matrix dimensions
#' @param boolean TRUE/FALSE. if true, matrix will only take values 0 and 1
#' @param slam if TRUE, will return matrix in simple_triplet_matrix form. If FALSE, will return SparseMatrix object
#' @return SparseMatrix or simple_triplet_matrix
#' @export
#' @examples
#'n.obs <- 1e5
#'n.vars <- 1e3
#'
#'#simulate a very sparse matrix (this matrix has many zeros and few ones)
#'x.sparse <- simSparseMatrix(sparsity = 0.99, dim = c(n.obs, n.vars), boolean = T)
simSparseMatrix <- function(sparsity, dim, boolean = T, slam = F) {
  #  generate a sparse X matrix with given sparsity
  #  and dimensions. can be boolean-valued or continuous
  stopifnot(length(dim) == 2)
  tot.cells <- dim[1] * dim[2]
  num.nonzero <- floor((1 - sparsity) * tot.cells)
  rows <- sample(1:dim[1], floor(num.nonzero + 0.01 * num.nonzero), replace = T)
  cols <- sample(1:dim[2], floor(num.nonzero + 0.01 * num.nonzero), replace = T)
  rc <- paste(rows, cols)
  nd.idx <- !duplicated(rc)
  if (sum(nd.idx) < num.nonzero) {
    rows <- sample(1:dim[1], floor(num.nonzero + 0.1 * num.nonzero), replace = T)
    cols <- sample(1:dim[2], floor(num.nonzero + 0.1 * num.nonzero), replace = T)
    rc <- paste(rows, cols)
    nd.idx <- !duplicated(rc)
  }
  if (sum(nd.idx) < num.nonzero) {
    rows <- sample(1:dim[1], floor(num.nonzero + 0.5 * num.nonzero), replace = T)
    cols <- sample(1:dim[2], floor(num.nonzero + 0.5 * num.nonzero), replace = T)
    rc <- paste(rows, cols)
    nd.idx <- !duplicated(rc)
  }
  rows <- rows[nd.idx][1:num.nonzero]
  cols <- cols[nd.idx][1:num.nonzero]
  if (boolean) {
    x.nonzero.vals <- 1
  } else {
    x.nonzero.vals <- rnorm(num.nonzero)
  }
  if (slam) {
    if (length(x.nonzero.vals) == 1) {x.nonzero.vals <- rep(x.nonzero.vals, num.nonzero)}
    sparse.mat <- simple_triplet_matrix(i = rows, j = cols, v = x.nonzero.vals, nrow = dim[1], ncol = dim[2])
  } else {
    sparse.mat <- sparseMatrix(i = rows, j = cols, x = x.nonzero.vals, dims = dim)
  }
  sparse.mat
}
