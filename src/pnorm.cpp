

#include "pnorm.h"


RcppExport SEXP ppnorm(SEXP x) 
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::Map;
    typedef Eigen::Map<VectorXd> MapVecd;
    
    
    const MapVecd xx(as<MapVecd>(x));
    const int n(xx.size());
    VectorXd ret(n);
    boost::math::normal_distribution<> std_normal(0.0, 1.0);
    
    for (int i = 0; i < n; i++) {
      ret(i) = cdf(std_normal, xx(i));
    }
    
    return wrap(ret);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}