
  
setGeneric("solveEigen", function(A, b, maxit = 500L, tol = 1e-5, method = c("BiCGSTAB", "CG")) {
  stopifnot(is.numeric(A) | inherits(A, "CsparseMatrix"))
  stopifnot(is.numeric(b))
  method <- match.arg(method)
  if (inherits(A, "CsparseMatrix")) {
    switch(method, 
           BiCGSTAB = .Call("BiCGSTAB_sparse_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions"),
           CG = .Call("conjugate_gradient_sparse_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions"))
  } else {
    stopifnot(is.matrix(A))
    switch(method,
           BiCGSTAB = .Call("BiCGSTAB_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions"),
           CG = .Call("conjugate_gradient_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions"))
  }
})

#' Solve the linear system Ax = b using conjugate gradient
#'
#' @param A square matrix input
#' @param b numeric vector. 
#' @param maxit integer. Maximum number of iterations
#' @param tol tolerance level
#' @return list
#' @export
#' @examples
#'n <- 2000
#'p <- 10
#'lambda <- 0.1
#'
#'set.seed(100)
#'x <- matrix(rnorm(n * p), ncol = p)
#'y <- x %*% rnorm(p) + rnorm(p)

#'b <- y
#'A <- tcrossprod(x) + lambda * diag(n)
#'
#'system.time(alpha.true <- solve(A, b))
#'system.time(alpha.cg <- solveCG(A, b))
#'
#'max(abs(alpha.cg$x - alpha.true))
setGeneric("solveCG", function(A, b, maxit = 500L, tol = 1e-5) {
  stopifnot(inherits(A, "matrix") | inherits(A, "CsparseMatrix"))
  stopifnot(is.numeric(b))
  n <- nrow(A)
  p <- ncol(A)
  bl <- length(b)
  stopifnot(n == p & p == bl)
  if (inherits(A, "CsparseMatrix")) {
    stop("not supported yet")
    #switch(method, 
    #       BiCGSTAB = .Call("BiCGSTAB_sparse_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions"),
    #       CG = .Call("conjugate_gradient_sparse_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions"))
  } else {
    stopifnot(is.matrix(A))
    CG = .Call("conjugate_gradient", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions")
    CG$x <- drop(CG$x)
    CG
  }
})

setMethod("solveEigen", signature(A = "matrix", b = "numeric", maxit = "numeric", 
                                 tol = "numeric", method = "character"),
          function(A, b, maxit = 50L, tol = 1e-5, method = c("BiCGSTAB", "CG")) {
            maxit <- as.integer(maxit)
            method <- match.arg(method)
            res <- switch(method,
                          BiCGSTAB = .Call("BiCGSTAB_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions"),
                          CG = .Call("conjugate_gradient_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions"))
            res
          })

setMethod("solveEigen", signature(A = "CsparseMatrix", b = "numeric", maxit = "numeric", 
                                tol = "numeric", method = "character"),
          function(A, b, maxit = 50L, tol = 1e-5, method = c("BiCGSTAB", "CG")) {
            maxit <- as.integer(maxit)
            method <- match.arg(method)
            res <- switch(method, 
                          BiCGSTAB = .Call("BiCGSTAB_sparse_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions"),
                          CG = .Call("conjugate_gradient_sparse_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions"))
            res
          })


setMethod("solveEigen", signature(A = "dgeMatrix", b = "numeric", maxit = "numeric", 
                                tol = "numeric", method = "character"),
          function(A, b, maxit = 50L, tol = 1e-5, method = c("BiCGSTAB", "CG")) {
            maxit <- as.integer(maxit)
            method <- match.arg(method)
            res <- switch(method,
                          BiCGSTAB = .Call("BiCGSTAB_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions"),
                          CG = .Call("conjugate_gradient_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions"))
            res
          })
