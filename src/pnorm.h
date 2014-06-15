#ifndef _rfunctions_PNORM_H
#define _rfunctions_PNORM_H

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector> 
#include <functional> 
#include <algorithm> 
#include <iostream>
#include <cmath>
#include <boost/math/distributions.hpp>


using namespace Rcpp;
using namespace RcppEigen;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp:wrap;

VectorXd pnorm(const VectorXd& x);


#endif