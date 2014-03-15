#ifndef _rfunctions_UTILS_H
#define _rfunctions_UTILS_H

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


RcppExport SEXP crossprodcpp(SEXP);

RcppExport SEXP xpwx(SEXP, SEXP);
                               
RcppExport SEXP subcpp(SEXP, SEXP);

RcppExport SEXP addcpp(SEXP, SEXP);

RcppExport SEXP subSparsecpp(SEXP, SEXP);

RcppExport SEXP addSparsecpp(SEXP, SEXP);

#endif