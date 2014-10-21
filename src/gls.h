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

RcppExport SEXP gls(SEXP, SEXP, SEXP, SEXP, SEXP);

RcppExport SEXP gls_half_cg(SEXP, SEXP, SEXP, SEXP, SEXP);

RcppExport SEXP gls_direct(SEXP, SEXP, SEXP);

#endif