
#' Solve the linear system (A'A + lambda * I)x = Ab using cgls
#'
#' @param A matrix or CsparseMatrix input
#' @param b numeric vector or numeric matrix (corresponding to multiple right-hand sides). 
#' @param lambda numeric positive scalar for ridge penalty 
#' @param maxit integer. Maximum number of iterations
#' @param tol tolerance level
#' @param init initialization of estimate
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
#'system.time( alpha.cgls <- cgls(x, b) )
#'
#'max(abs(alpha.cgls$x - alpha.true))
setGeneric("cgls", 
           function(A, b, 
                    lambda = 0, 
                    maxit = 50L, 
                    tol = 1e-5,
                    init = if (is.null(dim(b))){numeric(length(b))} 
                    else {array(0, dim = dim(b))} ) {
             #callNextMethod(A, b, lambda, maxit, tol, init)
             #callGeneric()
             standardGeneric("cgls")
})

setMethod("cgls", signature(A = "matrix"),
          function(A, b, lambda = 0, maxit = 50L, tol = 1e-5, 
                   init = if (is.null(dim(b))){numeric(ncol(A))} 
                          else {array(0, dim = c(ncol(A), ncol(b)))} ) {
            db <- dim(b)
            di <- dim(init)
            da <- dim(A)
            use.cgls <- TRUE
            if (!is.null(db) & !is.null(di)) {
              stopifnot(db[2] == di[2])
              stopifnot(da[2] == di[1])
              stopifnot(db[1] == da[1])
              if (db[2] > 1) {
                use.cgls <- FALSE
              }              
            } else if (!is.null(db) & is.null(di)) {
              stopifnot(da[2] == length(init))
              stopifnot(da[1] == db[1])
            } else if (is.null(db) & !is.null(di)) {
              stopifnot(di[1] == da[2])
              stopifnot(da[1] == length(b))
            } else {
              stopifnot(length(b) == da[1])
              stopifnot(da[2] == length(init))
            }

            if (use.cgls) {
              CG = .Call("cgls", 
                         A = A, 
                         b = drop(b), 
                         lambda = as.double(lambda),
                         maxit = as.integer(maxit), 
                         tol = as.double(tol),
                         init = drop(init), 
                         PACKAGE = "rfunctions")
            } else {
              CG = .Call("block_cgls", 
                         A = A, 
                         b = b, 
                         lambda = as.double(lambda),
                         maxit = as.integer(maxit), 
                         tol = as.double(tol),
                         init = init, 
                         PACKAGE = "rfunctions")
            }
            CG$x <- drop(CG$x)
            CG
          })



setMethod("cgls", signature(A = "dgeMatrix"),
          function(A, b, lambda = 0, maxit = 50L, tol = 1e-5, 
                   init = if (is.null(dim(b))){numeric(ncol(A))} 
                   else {array(0, dim = c(ncol(A), ncol(b)))} ) {
            db <- dim(b)
            di <- dim(init)
            da <- dim(A)
            use.cgls <- TRUE
            if (!is.null(db) & !is.null(di)) {
              stopifnot(db[2] == di[2])
              stopifnot(da[2] == di[1])
              stopifnot(db[1] == da[1])
              if (db[2] > 1) {
                use.cgls <- FALSE
              }              
            } else if (!is.null(db) & is.null(di)) {
              stopifnot(da[2] == length(init))
              stopifnot(da[1] == db[1])
            } else if (is.null(db) & !is.null(di)) {
              stopifnot(di[1] == da[2])
              stopifnot(da[1] == length(b))
            } else {
              stopifnot(length(b) == da[1])
              stopifnot(da[2] == length(init))
            }
            
            if (use.cgls) {
              CG = .Call("cgls", 
                         A = A, 
                         b = drop(b), 
                         lambda = as.double(lambda),
                         maxit = as.integer(maxit), 
                         tol = as.double(tol),
                         init = drop(init), 
                         PACKAGE = "rfunctions")
            } else {
              CG = .Call("block_cgls", 
                         A = A, 
                         b = b, 
                         lambda = as.double(lambda),
                         maxit = as.integer(maxit), 
                         tol = as.double(tol),
                         init = init, 
                         PACKAGE = "rfunctions")
            }
            CG$x <- drop(CG$x)
            CG
          })




setMethod("cgls", signature(A = "CsparseMatrix"),
          function(A, b, lambda = 0, maxit = 50L, tol = 1e-5, 
                   init = if (is.null(dim(b))){numeric(ncol(A))} 
                   else {array(0, dim = c(ncol(A), ncol(b)))} ) {
            db <- dim(b)
            di <- dim(init)
            da <- dim(A)
            use.cgls <- TRUE
            if (!is.null(db) & !is.null(di)) {
              stopifnot(db[2] == di[2])
              stopifnot(da[2] == di[1])
              stopifnot(db[1] == da[1])
              if (db[2] > 1) {
                use.cgls <- FALSE
              }              
            } else if (!is.null(db) & is.null(di)) {
              stopifnot(da[2] == length(init))
              stopifnot(da[1] == db[1])
            } else if (is.null(db) & !is.null(di)) {
              stopifnot(di[1] == da[2])
              stopifnot(da[1] == length(b))
            } else {
              stopifnot(length(b) == da[1])
              stopifnot(da[2] == length(init))
            }
            
            if (use.cgls) {
              CG = .Call("cgls_sparse", 
                         A = A, 
                         b = drop(b), 
                         lambda = as.double(lambda),
                         maxit = as.integer(maxit), 
                         tol = as.double(tol),
                         init = drop(init), 
                         PACKAGE = "rfunctions")
            } else {
              CG = .Call("block_cgls_sparse", 
                         A = A, 
                         b = b, 
                         lambda = as.double(lambda),
                         maxit = as.integer(maxit), 
                         tol = as.double(tol),
                         init = init, 
                         PACKAGE = "rfunctions")
            }
            CG$x <- drop(CG$x)
            CG
          })