
#' subtract b from a
#'
#' @param a matrix input
#' @param b matrix input
#' @return A matrix equal to b - a
#' @export
subtract <- function(a, b) a - b

#' add b to a
#'
#' @param a matrix input
#' @param b matrix input
#' @return A matrix equal to a + b
#' @export
add <- function(a, b) a + b

setMethod("subtract", signature(a = "dgCMatrix", b = "dgCMatrix"),
          function(a, b)
            .Call("subSparsecpp", BB = a, CC = b, PACKAGE = "rfunctions"),
          valueClass = "dgCMatrix")

setMethod("subtract", signature(a = "matrix", b = "matrix"),
          function(a, b)
            .Call("subcpp", BB = a, CC = b, PACKAGE = "rfunctions"),
          valueClass = "matrix")

setMethod("subtract", signature(a = "dsCMatrix", b = "dsCMatrix"),
          function(a, b)
            .Call("subSparsecpp", BB = as(a, "dgCMatrix"), CC = as(b, "dgCMatrix"), PACKAGE = "rfunctions"),
          valueClass = "dgCMatrix")

setMethod("subtract", signature(a = "dgCMatrix", b = "dsCMatrix"),
          function(a, b)
            .Call("subSparsecpp", BB = a, CC = as(b, "dgCMatrix"), PACKAGE = "rfunctions"),
          valueClass = "dgCMatrix")

setMethod("subtract", signature(a = "dsCMatrix", b = "dgCMatrix"),
          function(a, b)
            .Call("subSparsecpp", BB = as(a, "dgCMatrix"), CC = b, PACKAGE = "rfunctions"),
          valueClass = "dgCMatrix")

### ADDITION

setMethod("add", signature(a = "dgCMatrix", b = "dgCMatrix"),
          function(a, b)
            .Call("addSparsecpp", BB = a, CC = b, PACKAGE = "rfunctions"),
          valueClass = "dgCMatrix")


setMethod("add", signature(a = "matrix", b = "matrix"),
          function(a, b)
            .Call("addcpp", BB = a, CC = b, PACKAGE = "rfunctions"),
          valueClass = "matrix")

setMethod("add", signature(a = "dsCMatrix", b = "dsCMatrix"),
          function(a, b)
            .Call("addSparsecpp", BB = as(a, "dgCMatrix"), CC = as(b, "dgCMatrix"), PACKAGE = "rfunctions"),
          valueClass = "dgCMatrix")

setMethod("add", signature(a = "dgCMatrix", b = "dsCMatrix"),
          function(a, b)
            .Call("addSparsecpp", BB = a, CC = as(b, "dgCMatrix"), PACKAGE = "rfunctions"),
          valueClass = "dgCMatrix")

setMethod("add", signature(a = "dsCMatrix", b = "dgCMatrix"),
          function(a, b)
            .Call("addSparsecpp", BB = as(a, "dgCMatrix"), CC = b, PACKAGE = "rfunctions"),
          valueClass = "dgCMatrix")



setMethod("add", signature(a = "dgeMatrix", b = "dgeMatrix"),
          function(a, b)
            .Call("addcpp", BB = as(a, "matrix"), CC = as(b, "matrix"), PACKAGE = "rfunctions"),
          valueClass = "matrix")

setMethod("add", signature(a = "matrix", b = "dgeMatrix"),
          function(a, b)
            .Call("addcpp", BB = a, CC = as(b, "matrix"), PACKAGE = "rfunctions"),
          valueClass = "matrix")

setMethod("add", signature(a = "dgeMatrix", b = "matrix"),
          function(a, b)
            .Call("addcpp", BB = as(a, "matrix"), CC = b, PACKAGE = "rfunctions"),
          valueClass = "matrix")

setMethod("add", signature(a = "dgeMatrix", b = "numeric"),
          function(a, b)
            .Call("addcpp", BB = as(a, "matrix"), CC = b, PACKAGE = "rfunctions"),
          valueClass = "matrix")

setMethod("add", signature(a = "dgCMatrix", b = "numeric"),
          function(a, b)
            .Call("addcpp", BB = as(a, "matrix"), CC = b, PACKAGE = "rfunctions"),
          valueClass = "matrix")

setMethod("add", signature(a = "matrix", b = "dgCMatrix"),
          function(a, b)
            .Call("addcpp", BB = a, CC = as(b, "matrix"), PACKAGE = "rfunctions"),
          valueClass = "matrix")

setMethod("add", signature(a = "dgCMatrix", b = "matrix"),
          function(a, b)
            .Call("addcpp", BB = as(a, "matrix"), CC = b, PACKAGE = "rfunctions"),
          valueClass = "matrix")

setMethod("add", signature(a = "dgeMatrix", b = "dgCMatrix"),
          function(a, b)
            .Call("addcpp", BB = as(a, "matrix"), CC = as(b, "matrix"), PACKAGE = "rfunctions"),
          valueClass = "matrix")

setMethod("add", signature(a = "dgCMatrix", b = "dgeMatrix"),
          function(a, b)
            .Call("addcpp", BB = as(a, "matrix"), CC = as(b, "matrix"), PACKAGE = "rfunctions"),
          valueClass = "matrix")

setMethod("add", signature(a = "dgeMatrix", b = "dsCMatrix"),
          function(a, b)
            .Call("addcpp", BB = as(a, "matrix"), CC = as(b, "matrix"), PACKAGE = "rfunctions"),
          valueClass = "matrix")

setMethod("add", signature(a = "dsCMatrix", b = "dgeMatrix"),
          function(a, b)
            .Call("addcpp", BB = as(a, "matrix"), CC = as(b, "matrix"), PACKAGE = "rfunctions"),
          valueClass = "matrix")
