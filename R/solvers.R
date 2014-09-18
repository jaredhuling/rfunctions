
  
setGeneric("BiCGSTAB", function(A, b, maxit = 500L, tol = 1e-5) {
  stopifnot(is.numeric(A) | inherits(A, "CsparseMatrix"))
  stopifnot(is.numeric(b))
  if (inherits(x, "CsparseMatrix")) {
    .Call("BiCGSTAB_sparse_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions")
  } else {
    stopifnot(is.matrix(A))
    .Call("BiCGSTAB_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions")
  }
})

setMethod("BiCGSTAB", signature(A = "matrix", b = "numeric", maxit = "numeric", 
                                 tol = "numeric"),
          function(A, b, maxit = 50L, tol = 1e-5) {
            maxit <- as.integer(maxit)
            res <- .Call("BiCGSTAB_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions")
            res
          })

setMethod("BiCGSTAB", signature(A = "CsparseMatrix", b = "numeric", maxit = "numeric", 
                                tol = "numeric"),
          function(A, b, maxit = 50L, tol = 1e-5) {
            maxit <- as.integer(maxit)
            res <- .Call("BiCGSTAB_sparse_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions")
            res
          })


setMethod("BiCGSTAB", signature(A = "dgeMatrix", b = "numeric", maxit = "numeric", 
                                tol = "numeric"),
          function(A, b, maxit = 50L, tol = 1e-5) {
            maxit <- as.integer(maxit)
            res <- .Call("BiCGSTAB_eigen", A = A, b = b, maxit = maxit, tol = tol, PACKAGE = "rfunctions")
            res
          })
