

#' Solve the linear system A'Ax = Ab using lsmr
#'
#' @param A square matrix input
#' @param b numeric vector. 
#' @param maxit integer. Maximum number of iterations
#' @param tol tolerance level
#' @return list
#' @export
#' @examples
#'n <- 10000
#'p <- 100
#'
#'set.seed(100)
#'x <- simSparseMatrix(0.8, boolean = FALSE, dim = c(n, p)) 
#'b <- as.vector(x %*% rnorm(p) + rnorm(n))
#'
#'A <- crossprod(x)
#'xtb <- as.vector(crossprod(x, b))
#'
#'system.time( alpha.true <- solve(A, b) )
#'system.time( alpha.lsmr <- lsmr(x, b) )
#'
#'max(abs(alpha.lsmr$x - alpha.true))
setGeneric("lsmr", function(A, b, lambda = 0, atol = 1e-5, btol = 1e-5,
                            conlim = 1e18, maxit = 50L, localSize = 0) {
  stopifnot(inherits(A, "matrix") | inherits(A, "CsparseMatrix"))
  #stopifnot(is.numeric(b))
  n <- nrow(A)
  p <- ncol(A)
  bl <- length(b)
  if (inherits(A, "CsparseMatrix")) {
    res = .Call("lsmr_sparse", A = A, b = b, lambda = lambda, 
                atol = atol, btol = btol, 
                conlim = conlim, maxit = maxit,
                localSize = localSize, PACKAGE = "rfunctions")
    res$x <- drop(res$x)
    res    
  } else {
    stopifnot(is.matrix(A))
    stop("not implemented yet")
  }
})

