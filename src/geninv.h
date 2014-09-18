#ifndef _rfunctions_GENINV_H
#define _rfunctions_GENINV_H

//#define EIGEN_SUPERLU_SUPPORT
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/LU>
#include <Eigen/SparseCholesky>
//#include <Eigen/SuperLUSupport>
//#include <Eigen/IterativeSolvers>
#include <vector> 
#include <functional> 
#include <algorithm> 
#include <iostream>
#include <cmath>



RcppExport SEXP geninv(SEXP);

//RcppExport SEXP geninv_sparse(SEXP);

RcppExport SEXP ginv(SEXP);

#endif