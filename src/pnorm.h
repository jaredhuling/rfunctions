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
#include <boost/random/normal_distribution.hpp> 
#include <boost/math/distributions.hpp>

using Eigen::MatrixXd;
using Eigen::VectorXd;

VectorXd pnorm(const VectorXd& x) {



#endif