

big.crossprod <- function(x)
{
  if (!is.big.matrix(x))
  {
    stop("object must be a big.matrix object")
  }
  .Call("crossprod_big", x@address)
}

big.colSums <- function(x)
{
  if (!is.big.matrix(x))
  {
    stop("object must be a big.matrix object")
  }
  .Call("colsums_big", x@address)
}