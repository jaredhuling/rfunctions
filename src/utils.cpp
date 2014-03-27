

#include "utils.h"


//computes X'WX where W is diagonal (input w as vector)
MatrixXd xtwx(const MatrixXd& xx, const MatrixXd& ww) {
  const int n(xx.cols());
  MatrixXd AtWA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx.adjoint() * ww.asDiagonal()));
  return (AtWA);
}

MatrixXd xtx(const MatrixXd& xx) {
  const int n(xx.cols());
  MatrixXd AtA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx.adjoint()));
  return (AtA);
}

MatrixXd xxt(const MatrixXd& xx) {
  const int m(xx.rows());
  MatrixXd AtA(MatrixXd(m, m).setZero().
    selfadjointView<Lower>().rankUpdate(xx));
  return (AtA);
}