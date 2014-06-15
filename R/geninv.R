

#' Compute The Moore-Penrose generalized inverse of a matrix
#'
#' @param A matrix 
#' @return The pseudoinverse of matrix A.
#' @references Courrieu (2005) Fast Computation of Moore-Penrose Inverse Matrices, http://arxiv.org/ftp/arxiv/papers/0804/0804.4809.pdf 
#' Neural Information Processing - Letters and Reviews, Vol. 8(2), Aug. 2005
#' @export
#' @examples
#'n <- 1000
#'p <- 500
#'
#'x <- matrix(rnorm(n * (p - 1)), n, p-1)
#'x <- cbind(x, rowMeans(x))
#'
#'## compute X'X
#'xpx <- crossprodcpp(x)
#'
#'## compute generalized inverse of X'X
#'inv <- geninv(xpx)
#'
#'## check if we have computed a generalized inverse
#'all.equal(xpx, xpx %*% inv %*% xpx)
setGeneric("geninv", function(A) {
  stopifnot(is.numeric(A) | inherits(A, "CsparseMatrix"))

  if (inherits(A, "CsparseMatrix")) {
    .Call("geninv_sparse", GG = A, PACKAGE = "rfunctions")
  } else {
    .Call("geninv", GG = A, PACKAGE = "rfunctions")
  }
})