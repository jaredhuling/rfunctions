#ifndef _rfunctions_LSMR_H
#define _rfunctions_LSMR_H

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

RcppExport SEXP lsmr_sparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


#endif