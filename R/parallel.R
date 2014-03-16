

runParallel <- function(ncores = NULL){
  require(parallel)
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