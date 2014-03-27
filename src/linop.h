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

inline Eigen::ArrayXd Dplus(const Eigen::ArrayXd& d);

RcppExport SEXP crossprodcpp(SEXP);

RcppExport SEXP xpwx(SEXP, SEXP);
                               
RcppExport SEXP subcpp(SEXP, SEXP);

RcppExport SEXP addcpp(SEXP, SEXP);

RcppExport SEXP subSparsecpp(SEXP, SEXP);

RcppExport SEXP addSparsecpp(SEXP, SEXP);

RcppExport SEXP fastRank(SEXP);

#endif