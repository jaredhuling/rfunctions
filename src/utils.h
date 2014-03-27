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

//computes X'X 
MatrixXd xtx(const MatrixXd& xx);

//computes XX' 
MatrixXd xxt(const MatrixXd& xx);



#endif