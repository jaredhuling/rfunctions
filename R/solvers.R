
  
setGeneric("solveEigen", function(A, b, maxit = 500L, tol = 1e-5, method = c("BiCGSTAB", "CG")) {
  stopifnot(is.numeric(A) | inherits(A, "CsparseMatrix"))
  stopifnot(is.numeric(b))
  method <- match.arg(method)
  if (inherits(x, "CsparseMatrix")) {
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
