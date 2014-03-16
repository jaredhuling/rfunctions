

#' Compute Largest Singular Value of Matrix x
#'
#' @param x matrix input
#' @param v numeric vector. Initialize GKL bidiagonalization with random vector on unit sphere
#' @param maxit integer. Maximum number of Lanczos iterations
#' @param reorthog. Takes values in 0:2. 0 for no reorthogonalization, 1 for reorthogonalization of 
#' V vectors (slower, more accurate), 2 for reorthogonalization of V and U vectors (slowest and most 
#' memory, most accurate. Not implemented yet)
#' @param upper.bound.prob upper bound probability for the largest singular value
#' @return list 
#' @references Hochstenbach (2013) Probabilistic upper bounds for the matrix two-norm, http://link.springer.com/article/10.1007/s10915-013-9716-x 
#' Journal of Scientific Computing, Vol. 57(3), Dec. 2013
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
setGeneric("gklBidiag", function(x, v, maxit, reorthog = 0, upper.bound.prob = NULL) {
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(v))
  reorthog <- as.integer(reorthog)
  .Call("GKLBidiag", A = x, v = v, k = maxit, reorthog = reorthog, PACKAGE = "rfunctions")
})

setMethod("gklBidiag", signature(x = "matrix", v = "numeric", maxit = "numeric", 
                                 reorthog = "numeric", upper.bound.prob = "numeric"),
            function(x, v, maxit = 50L, reorthog = 0, upper.bound.prob = NULL) {
              maxit <- as.integer(maxit)
              reorthog <- as.integer(reorthog)
              delt <- computeDelta(eps = 1 - upper.bound.prob, p = ncol(x))
              res <- .Call("GKLBidiag", A = x, v = v, k = maxit, reorthog = reorthog, PACKAGE = "rfunctions")
              upper <- computeUpperBound(res, delt)
              res$upper.bound <- upper
              res$upper.bound.prob <- upper.bound.prob
              res
            })

setMethod("gklBidiag", signature(x = "dgeMatrix", v = "numeric", maxit = "numeric", 
                                 reorthog = "numeric", upper.bound.prob = "numeric"),
          function(x, v, maxit = 50L, reorthog = 0, upper.bound.prob = NULL) {
            maxit <- as.integer(maxit)
            reorthog <- as.integer(reorthog)
            delt <- computeDelta(eps = 1 - upper.bound.prob, p = ncol(x))
            res <- .Call("GKLBidiag", A = x, v = v, k = maxit, reorthog = reorthog, PACKAGE = "rfunctions")
            upper <- computeUpperBound(res, delt)
            res$upper.bound <- upper
            res$upper.bound.prob <- upper.bound.prob
            res
          })


setMethod("gklBidiag", signature(x = "CsparseMatrix", v = "numeric", maxit = "numeric", 
                                 reorthog = "numeric", upper.bound.prob = "numeric"),
          function(x, v, maxit = 50L, reorthog = 0, upper.bound.prob = NULL) {
            maxit <- as.integer(maxit)
            reorthog = as.integer(reorthog)
            delt <- computeDelta(eps = 1 - upper.bound.prob, p = ncol(x))
            res <- .Call("GKLBidiagSparse", A = x, v = v, k = maxit, reorthog = reorthog, PACKAGE = "rfunctions")
            upper <- computeUpperBound(res, delt)
            res$upper.bound <- upper
            res$upper.bound.prob <- upper.bound.prob
            res
          })


setMethod("gklBidiag", signature(x = "matrix", v = "numeric", maxit = "numeric", 
                                 reorthog = "numeric", upper.bound.prob = "missing"),
          function(x, v, maxit = 50L, reorthog = 0, upper.bound.prob = NULL) {
            maxit <- as.integer(maxit)
            reorthog <- as.integer(reorthog)
            res <- .Call("GKLBidiag", A = x, v = v, k = maxit, reorthog = reorthog, PACKAGE = "rfunctions")
            res$upper.bound <- NULL
            res$upper.bound.prob <- NULL
            res
          })

setMethod("gklBidiag", signature(x = "dgeMatrix", v = "numeric", maxit = "numeric", 
                                 reorthog = "numeric", upper.bound.prob = "missing"),
          function(x, v, maxit = 50L, reorthog = 0, upper.bound.prob = NULL) {
            maxit <- as.integer(maxit)
            reorthog <- as.integer(reorthog)
            res <- .Call("GKLBidiag", A = x, v = v, k = maxit, reorthog = reorthog, PACKAGE = "rfunctions")
            res$upper.bound <- NULL
            res$upper.bound.prob <- NULL
            res
          })


setMethod("gklBidiag", signature(x = "CsparseMatrix", v = "numeric", maxit = "numeric", 
                                 reorthog = "numeric", upper.bound.prob = "missing"),
          function(x, v, maxit = 50L, reorthog = 0, upper.bound.prob = NULL) {
            maxit <- as.integer(maxit)
            reorthog <- as.integer(reorthog)
            res <- .Call("GKLBidiagSparse", A = x, v = v, k = maxit, reorthog = reorthog, PACKAGE = "rfunctions")
            res$upper.bound <- NULL
            res$upper.bound.prob <- NULL
            res
          })

