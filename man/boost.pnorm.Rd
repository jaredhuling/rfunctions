\name{boost.pnorm}
\alias{boost.pnorm}
\title{distribution function}
\usage{
boost.pnorm(q, mean = 0, sd = 1, lower.tail = TRUE)
}
\arguments{
  \item{q}{vector of quantiles}

  \item{mean}{scalar mean (vector will be implemented)}

  \item{sd}{scalar standard deviation (vector will be
  implemented)}

  \item{lower.tail}{logical; if TRUE (default),
  probabilities are P[X <= x] otherwise, P[X > x].}
}
\value{
the distribution function
}
\description{
distribution function
}
\examples{
library(rfunctions)
vals <- rnorm(100, mean = 4, sd = 2)
boost.pnorm(vals, mean = 4, sd = 2, lower.tail = FALSE)
}

