

#' Compute Largest Singular Value of Matrix x
#'
#' @param x matrix input
#' @param v numeric vector. Initialize GKL bidiagonalization with random vector on unit sphere
#' @param maxit integer. Maximum number of Lanczos iterations
#' @param upper.bound.prob upper bound probability for the largest singular value
#' @return list 
#' @export
#' @examples
#'n.obs <- 1e5
#'n.vars <- 150
#'
#'x <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#'
#'## compute largest singular value of x
#'lanczos <- gklBidiag(x, runif(ncol(x)), 10L, 0.99)
#'str(lanczos)
setGeneric("gklBidiag", function(x, v, maxit, upper.bound.prob = NULL) {
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(v))
  .Call("GKLBidiag", A = x, v = v, k = maxit, PACKAGE = "rfunctions")
})

setMethod("gklBidiag", signature(x = "matrix", v = "numeric", maxit = "numeric", upper.bound.prob = "numeric"),
            function(x, v, maxit = 50L, upper.bound.prob = NULL) {
              maxit <- as.integer(maxit)
              delt <- computeDelta(eps = 1 - upper.bound.prob, p = ncol(x))
              res <- .Call("GKLBidiag", A = x, v = v, k = maxit, PACKAGE = "rfunctions")
              upper <- computeUpperBound(res, delt)
              res$upper.bound <- upper
              res$upper.bound.prob <- upper.bound.prob
              res
            })

setMethod("gklBidiag", signature(x = "dgeMatrix", v = "numeric", maxit = "numeric", upper.bound.prob = "numeric"),
          function(x, v, maxit = 50L, upper.bound.prob = NULL) {
            maxit <- as.integer(maxit)
            delt <- computeDelta(eps = 1 - upper.bound.prob, p = ncol(x))
            res <- .Call("GKLBidiag", A = x, v = v, k = maxit, PACKAGE = "rfunctions")
            upper <- computeUpperBound(res, delt)
            res$upper.bound <- upper
            res$upper.bound.prob <- upper.bound.prob
            res
          })


setMethod("gklBidiag", signature(x = "CsparseMatrix", v = "numeric", maxit = "numeric", upper.bound.prob = "numeric"),
          function(x, v, maxit = 50L, upper.bound.prob = NULL) {
            maxit <- as.integer(maxit)
            delt <- computeDelta(eps = 1 - upper.bound.prob, p = ncol(x))
            res <- .Call("GKLBidiagSparse", A = x, v = v, k = maxit, PACKAGE = "rfunctions")
            upper <- computeUpperBound(res, delt)
            res$upper.bound <- upper
            res$upper.bound.prob <- upper.bound.prob
            res
          })


setMethod("gklBidiag", signature(x = "matrix", v = "numeric", maxit = "numeric", upper.bound.prob = "missing"),
          function(x, v, maxit = 50L, upper.bound.prob = NULL) {
            maxit <- as.integer(maxit)
            res <- .Call("GKLBidiag", A = x, v = v, k = maxit, PACKAGE = "rfunctions")
            res$upper.bound <- NULL
            res$upper.bound.prob <- NULL
            res
          })

setMethod("gklBidiag", signature(x = "dgeMatrix", v = "numeric", maxit = "numeric", upper.bound.prob = "missing"),
          function(x, v, maxit = 50L, upper.bound.prob = NULL) {
            maxit <- as.integer(maxit)
            res <- .Call("GKLBidiag", A = x, v = v, k = maxit, PACKAGE = "rfunctions")
            res$upper.bound <- NULL
            res$upper.bound.prob <- NULL
            res
          })


setMethod("gklBidiag", signature(x = "CsparseMatrix", v = "numeric", maxit = "numeric", upper.bound.prob = "missing"),
          function(x, v, maxit = 50L, upper.bound.prob = NULL) {
            maxit <- as.integer(maxit)
            res <- .Call("GKLBidiagSparse", A = x, v = v, k = maxit, PACKAGE = "rfunctions")
            res$upper.bound <- NULL
            res$upper.bound.prob <- NULL
            res
          })

