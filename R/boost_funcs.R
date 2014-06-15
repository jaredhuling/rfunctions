
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