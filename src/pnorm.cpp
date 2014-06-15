

#include "pnorm.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;

VectorXd pnorm(const VectorXd& x) {
  const int n(xx.size());
  VectorXd ret(VectorXd(n));
  boost::math::normal_distribution<> std_normal(0.0, 1.0);
  
  for (int i = 0, i < n, i++) {
    ret(i) = cdf(std_normal, x(i));
  }
  return (ret);
}