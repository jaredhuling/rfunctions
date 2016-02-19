

#' Compute X'X or X'WX where W is diagonal weight matrix
#'
#' @param x matrix input
#' @param w numeric vector. If specified, X'diag(w)X will be computed
#' @return A symmetric matrix X'X
#' @export
#' @examples
#'n.obs <- 1e5
#'n.vars <- 150
#'
#'x <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#'
#'## compute X'X
#'xpx <- crossprodcpp(x)
#'
#'weights <- runif(n.obs)
#'
#'## compute X'WX
#'xpwx <- crossprodcpp(x, weights)
setGeneric("crossprodcpp", function(x, w = NULL) crossprod(x))

#crossprodcpp <- function(x, w = NULL) {
#  .Call("crossprodcpp", X = x, PACKAGE = "rfunctions")
#}

setMethod("crossprodcpp", signature(x = "dgeMatrix", w = "missing"),
          function(x) .Call("crossprodcpp", X = x, PACKAGE = "rfunctions"),
          valueClass = "dgeMatrix")

setMethod("crossprodcpp", signature(x = "matrix", w = "missing"),
          function(x) .Call("crossprodcpp", X = x, PACKAGE = "rfunctions"),
          valueClass = "matrix")

setMethod("crossprodcpp", signature(x = "dgCMatrix", w = "missing"),
          function(x) forceSymmetric(.Call("crossprodSparsecpp", X = x, PACKAGE = "rfunctions"), uplo="U"),
          valueClass = "dscMatrix")

setMethod("crossprodcpp", signature(x = "dsCMatrix", w = "missing"),
          function(x) forceSymmetric(.Call("crossprodSparsecpp", X = x, PACKAGE = "rfunctions"), uplo="U"),
          valueClass = "dscMatrix")

setMethod("crossprodcpp", signature(x = "dgeMatrix", w = "numeric"),
          function(x, w = NULL) .Call("xpwx", X = x, W = w, PACKAGE = "rfunctions"),
          valueClass = "dgeMatrix")

setMethod("crossprodcpp", signature(x = "matrix", w = "numeric"),
          function(x, w = NULL) .Call("xpwx", X = x, W = w, PACKAGE = "rfunctions"),
          valueClass = "matrix")


#' Compute X'X in chunks
#'
#' @param X matrix input
#' @param row.chunk integer. Number of rows per chunk. Specify NULL for no chunking
#' @param sparse Boolean value to specify if matrix X is sparse 
#' @return A symmetric matrix X'X
#' @export
#' @examples
#' crossprodChunk(matrix(rnorm(1e4 * 1e2), 1e4, 1e2), row.chunk = 500, FALSE)
crossprodChunk <- function(X, row.chunk = NULL, sparse = FALSE){
  if (sparse) {
    if (!is.null(w)) X<-sqrt(w) * X
    if (class(X) != "dgCMatrix")
      X <- as(X, "dgCMatrix")
    new.B <- Matrix:::crossprod(X)
  }
  else {
    if (is.null(row.chunk)) {
      new.B <- .Call("crossprodcpp", X = X, PACKAGE = "rfunctions")
    } else {
      if (row.chunk >= (nrow(X) - 1)) {
        row.chunk <- nrow(X)
      }
      if (row.chunk <= 1) {
        row.chunk <- 2
      }
      mod <- nrow(X)%%row.chunk
      last.block <- (mod > 0)
      G <- nrow(X)%/%row.chunk - (last.block * (2 * row.chunk <=
                                                  nrow(X)))
      new.B <- matrix(0, ncol(X), ncol(X))
      a <- row.chunk
      j <- 1
      for (i in 1:G) {
        B <- .Call("crossprodcpp", X = X[j:(i * a), ], PACKAGE = "rfunctions") # crossprod(X[j:(i * a), ]) #
        new.B <- new.B + B
        j <- j + a
      }
      if (mod > 0) {
        new.B <- new.B + .Call("crossprodcpp", X = X[j:nrow(X), ], PACKAGE = "rfunctions") # crossprod(X[j:nrow(X), ]) #
      }
    }
  }
  as(new.B, "matrix")
}


#call fast matrix rank function
fastRank <- function(M) {
  stopifnot(inherits(M, "matrix"))
  stopifnot(is.numeric(M))
  .Call("fastRank", AA = M, PACKAGE = "rfunctions")
}


#' Compute X[,indices]'y 
#'
#' @param x matrix input
#' @param y numeric vector. 
#' @param indices integer vector of X columns to use
#' @return x[,indices]'y
#' @export
#' @examples
#'n.obs <- 1e4
#'n.vars <- 2000
#'n.idx <- 500
#'
#'x <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#'y <- rnorm(n.obs)
#'idx <- sample.int(n.vars, n.idx)
#'
#'## compute x[,idx]'y
#'xpy <- sliced.crossprod(x, y, idx)
#'
#'xpy2 <- crossprod(x[,idx], y)
#'
#'all.equal(xpy, xpy2)
#'
#'@export
sliced.crossprod <- function(x, y, indices) {
  stopifnot(inherits(x, "matrix"))
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(y))
  .Call("sliced_crossprod", X_ = x, Y_ = y, idx_ = as.integer(indices), PACKAGE = "rfunctions")
}

