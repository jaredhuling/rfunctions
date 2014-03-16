\name{runParallel}
\alias{runParallel}
\title{Initialized parallel workers}
\usage{
runParallel(ncores = 1L)
}
\arguments{
  \item{ncores}{integer number of cores to initialize}
}
\value{
cluster object
}
\description{
Initialized parallel workers
}
\examples{
cl <- runParallel(ncores = 4L)
x <- 0
foreach(i=1:100) \%dopar$ {
x <- x + i ^ 2
}
}

