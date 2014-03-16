

#ifndef _rfunctions_LANCZOS_H
#define _rfunctions_LANCZOS_H

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


RcppExport SEXP GKLBidiag(SEXP, SEXP, SEXP);

RcppExport SEXP GKLBidiagSparse(SEXP, SEXP, SEXP);

RcppExport SEXP BidiagPoly(SEXP, SEXP, SEXP);

#endif