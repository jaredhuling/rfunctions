#ifndef _rfunctions_LINOP_H
#define _rfunctions_LINOP_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/SVD>
#include <vector> 
#include <functional> 
#include <algorithm> 
#include <iostream>
#include <cmath>

using namespace Rcpp;
using namespace RcppEigen;


RcppExport SEXP crossprod_big(SEXP);

RcppExport SEXP colsums_big(SEXP);

#endif