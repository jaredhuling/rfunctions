#ifndef _rfunctions_TESTING_H
#define _rfunctions_TESTING_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <igl/slice.h>
#include <Eigen/SVD>
#include <vector> 
#include <functional> 
#include <algorithm> 
#include <iostream>
#include <cmath>
#include "Spectra/SymEigsSolver.h"
#include "MatOp.h"

using namespace Rcpp;
using namespace RcppEigen;

RcppExport SEXP crossprodeig(SEXP);

RcppExport SEXP matcopytest(SEXP);

RcppExport SEXP listtest(SEXP, SEXP, SEXP);

RcppExport SEXP tcrossprodvec (SEXP);

RcppExport SEXP crossprodint(SEXP);

RcppExport SEXP crossprodrowmajor(SEXP, SEXP);

RcppExport SEXP lentest(SEXP);

RcppExport SEXP crossprodxval(SEXP, SEXP);

RcppExport SEXP dev_resid_logistic(SEXP, SEXP);

RcppExport SEXP crossprodsubset(SEXP, SEXP);

RcppExport SEXP matveccrossprodidx(SEXP, SEXP, SEXP);

RcppExport SEXP matvecprodidx(SEXP, SEXP, SEXP);

RcppExport SEXP appendMat(SEXP, SEXP);

#endif