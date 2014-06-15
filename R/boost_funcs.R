

boost.pnorm <- function(q, mean = 0, sd = 1, lower.tail = TRUE) {
  .Call("ppnorm", q = q, mean = mean, sd = sd, lower.tail = lower.tail)
}