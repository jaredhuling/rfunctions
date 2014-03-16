

cosf <- function(t, p) {
  (cos(t)) ^ (p - 2)
}

zeroDelt <- function(t, eps, pp) {
  stopifnot(t <= 1)
  given <- eps * integrate(cosf, 0, pi / 2, p = pp)$value
  integrate(cosf, 0, asin(t), p = pp)$value - given
}

computeDelta <- function(eps, p) {
  #compute delta value for Singular Value upper bound. Ref: Hochstenbach (2013)
  1 / uniroot(zeroDelt, interval = c(0, 0.5), eps = eps, pp = p)$root
}

polyfunc <- function(t, alpha, beta, delta) {
  .Call("BidiagPoly", X = t^2, alpha = alpha, beta = beta) * t - delta
}
vpoly <- Vectorize(polyfunc, vectorize.args = "t")

computeUpperBound <- function(lanc.obj, delta) {
  uniroot(vpoly, lower = lanc.obj$d, upper = lanc.obj$d * 100, alpha = lanc.obj$alpha, 
          beta = lanc.obj$beta, delta = delta)$root
}

