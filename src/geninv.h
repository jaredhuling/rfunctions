#ifndef _rfunctions_GENINV_H
#define _rfunctions_GENINV_H


#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/LU>
#include <vector> 
#include <functional> 
#include <algorithm> 
#include <iostream>
#include <cmath>



RcppExport SEXP geninv(SEXP);



#endif