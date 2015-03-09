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

RcppExport SEXP BiCGSTAB_eigen(SEXP, SEXP, SEXP, SEXP);

RcppExport SEXP BiCGSTAB_sparse_eigen(SEXP, SEXP, SEXP, SEXP);

RcppExport SEXP conjugate_gradient_eigen(SEXP, SEXP, SEXP, SEXP);

RcppExport SEXP conjugate_gradient_sparse_eigen(SEXP, SEXP, SEXP, SEXP);

RcppExport SEXP conjugate_gradient(SEXP, SEXP, SEXP, SEXP);

RcppExport SEXP ILU_prec_conjugate_gradient_sparse(SEXP, SEXP, SEXP, SEXP);

RcppExport SEXP conjugate_gradient_sparse(SEXP, SEXP, SEXP, SEXP);

RcppExport SEXP block_conjugate_gradient(SEXP, SEXP, SEXP, SEXP);

RcppExport SEXP block_conjugate_gradient_sparse(SEXP, SEXP, SEXP, SEXP);

RcppExport SEXP incompleteLUT_eig(SEXP, SEXP);

#endif