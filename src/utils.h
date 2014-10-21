#ifndef _rfunctions_UTILS_H
#define _rfunctions_UTILS_H


#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector> 
#include <functional> 
#include <algorithm> 
#include <iostream>
#include <cmath>


using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::Lower;



//computes X'WX where W is diagonal (input w as vector)
MatrixXd xtwx(const MatrixXd& xx, const MatrixXd& ww);

//computes X'SX where S is not diagonal (input ss as matrix)
MatrixXd xtsx(const MatrixXd& xx, const MatrixXd& ss);

//computes X'X 
MatrixXd xtx(const MatrixXd& xx);

//computes XX' 
MatrixXd xxt(const MatrixXd& xx);

//solve Ax = b using conjugate gradient
MatrixXd conjugate_gradient(const MatrixXd& A, const VectorXd& b, int maxit, double tol);

#endif