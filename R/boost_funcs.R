
#'  distribution function
#'
#' @param q vector of quantiles
#' @param mean scalar mean (vector will be implemented)
#' @param sd scalar standard deviation (vector will be implemented)
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x] otherwise, P[X > x].
#' @return the distribution function
#' @export
#' @examples
#' library(rfunctions)
#' vals <- rnorm(100, mean = 4, sd = 2)
#' boost.pnorm(vals, mean = 4, sd = 2, lower.tail = FALSE)
boost.pnorm <- function(q, mean = 0, sd = 1, lower.tail = TRUE) {
  if (!(inherits(q, "numeric") | inherits(q, "integer"))) stop("Non-numeric argument to mathematical function")
  .Call("ppnorm", q = as.numeric(q), mean = as.numeric(mean), 
        sd = as.numeric(sd), lower.tail = lower.tail)
}


#'  distribution function
#'
#' @param n integer number of samples
#' @param mean vector mean of length equal to ncol(sigma) = nrow(sigma). If scalar, assumed the same for each dimension
#' @param sigma 
#' @param seed integer; if NULL (default), seed is not forced.
#' @return matrix of samples from multivariate normal distribution 
#' @export
#' @examples
#' library(rfunctions)
#' vals <- rnorm(100, mean = 4, sd = 2)
#' boost.pnorm(vals, mean = 4, sd = 2, lower.tail = FALSE)
multivarNorm <- function(n = 1L, mean = 0, sigma, seed = NULL) {
  p <- length(mean)
  ds <- dim(sigma)
  if (p == 1 & p < ds[1]) {
    mean <- rep(mean, ds[1])
  } else {
    if (!all(ds == c(p, p))) 
      stop("incompatible arguments")
  }
  
  
  if (!(inherits(n, "numeric") | inherits(n, "integer"))) stop("Non-numeric argument to mathematical function")
  n <- ceiling(n)
  stopifnot(is.numeric(sigma))
  
  if (is.null(seed)) {
    seed <- -999999
  }
  
  .Call("multivarNorm", 
        n = ceiling(n), 
        mean = as.numeric(mean), 
        sigma = sigma, 
        seed = seed)
}

