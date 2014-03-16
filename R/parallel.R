
#' Initialized parallel workers
#'
#' @param ncores integer number of cores to initialize
#' @return cluster object 
#' @export
#' @examples
#'cl <- runParallel(ncores = 4L)
#'x <- 0
#'foreach(i=1:100) %dopar$ {
#'x <- x + i ^ 2
#'}
runParallel <- function(ncores = 1L){
  require(parallel)
  ncores <- as.integer(ncores)
  switch(Sys.info()[["sysname"]], 
         "Windows" = {
           library(doSNOW)
           #library(doParallel)
           num.cores <- ifelse(is.null(ncores), detectCores() - 1, ncores)
           cl <- makePSOCKcluster(max(1, num.cores)) #; registerDoSNOW(cl)
           clusterCall(cl, worker.init)
           registerDoParallel(cl)
         },
         "Linux" = {
           library(doMC); registerDoMC()
           if(is.null(ncores)) {
             options(cores = detectCores() - 1)
           } else {options(cores = ncores)}
           cl <- NULL
         })
  cl
}