

#include "utils.h"


//computes X'WX where W is diagonal (input w as vector)
MatrixXd xtwx(const MatrixXd& xx, const MatrixXd& ww) {
  const int n(xx.cols());
  MatrixXd AtWA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx.adjoint() * ww.asDiagonal()));
  return (AtWA);
}

//computes X'SX where S is not diagonal (input ss as matrix)
MatrixXd xtsx(const MatrixXd& xx, const MatrixXd& ss) {
  const int n(xx.cols());
  MatrixXd AtSA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx.adjoint() * ss));
  return (AtSA);
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


MatrixXd conjugate_gradient(const MatrixXd& A, const VectorXd& b, int maxit, double tol)
{

  const int n(A.cols());
  VectorXd x(n);
  VectorXd r(n);
  VectorXd p(n);
  VectorXd Ap(n);
  x.fill(0);
  
  double rsold;
  double rsnew;
  double alpha;
  int iters = maxit; 
  
  r = b; 
  p = r;
  rsold = r.squaredNorm();
  
  for (int i = 0; i < maxit; i++) {
    Ap = A * p;
    alpha = rsold / (p.transpose() * Ap);
    x = x + (alpha * p.array()).matrix();
    r = r - (alpha * Ap.array()).matrix();
    rsnew = r.squaredNorm();
    if (sqrt(rsnew) < tol) {
      iters = i;
      break;
    }
    p = r + ((rsnew / rsold) * p.array()).matrix();
    rsold = rsnew;
  }
  return(x);
}


