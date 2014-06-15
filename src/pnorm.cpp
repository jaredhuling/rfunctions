

#include "pnorm.h"


RcppExport SEXP ppnorm(SEXP q, SEXP mean, SEXP sd, SEXP lower_tail) 
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::Map;
    using Rcpp::as;
    typedef Eigen::Map<VectorXd> MapVecd;
    
    
    const MapVecd qq(as<MapVecd>(q));
    const int n(qq.size());
    const double mean_(as<double>(mean));
    const double sd_(as<double>(sd));
    const bool lt(as<bool>(lower_tail));
    
    VectorXd ret(n);
    
    boost::math::normal_distribution<> std_normal(mean_, sd_);
    
    if (lt) {
      // lower tail probability
      for (int i = 0; i < n; i++) {
        ret(i) = cdf(std_normal, qq(i));
      }
    } else {
      // upper tail probability
      for (int i = 0; i < n; i++) {
        ret(i) = cdf(complement (std_normal, qq(i)) );
      }
    }
    
    return wrap(ret);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}