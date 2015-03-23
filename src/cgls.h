#ifndef _rfunctions_CGLS_H
#define _rfunctions_CGLS_H

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


RcppExport SEXP cgls(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

RcppExport SEXP block_cgls(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

RcppExport SEXP cgls_sparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

RcppExport SEXP block_cgls_sparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

#endif