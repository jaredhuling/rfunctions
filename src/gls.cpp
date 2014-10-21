#include <iostream>
#include "linop.h"
#include "utils.h"
#include "math.h"
#include <Eigen/Cholesky>

using namespace Rcpp;
using namespace RcppEigen;

RcppExport SEXP gls(SEXP X, SEXP S, SEXP Y, SEXP maxit, SEXP tol)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Rcpp::List;
    
    typedef Map<VectorXd> MapVecd;
    typedef Map<Eigen::MatrixXd> MapMatd;
    
    const int maxiter(as<int>(maxit));
    const double toler(as<double>(tol));
    const Eigen::Map<MatrixXd> XX(as<MapMatd>(X));
    const Eigen::Map<MatrixXd> SS(as<MapMatd>(S));
    const Eigen::Map<VectorXd> YY(as<MapVecd>(Y));
    
    const int n(XX.rows());
    const int p(XX.cols());
    
    MatrixXd SX(MatrixXd(n, p));
    //MatrixXd XSX(MatrixXd(p, p));
    
    for (int i = 0; i < p; i++) {
      SX.col(i) = conjugate_gradient(SS, XX.col(i), maxiter, toler);
    }
    
    MatrixXd XSX((XX.adjoint()) * SX);
    VectorXd XSY((SX.adjoint()) * YY);
    
    VectorXd beta = conjugate_gradient(XSX, XSY, maxiter, toler);
    
    
    return List::create(Named("beta") = beta);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

RcppExport SEXP gls_half_cg(SEXP X, SEXP S, SEXP Y, SEXP maxit, SEXP tol)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    
    using Rcpp::List;
    
    typedef Map<VectorXd> MapVecd;
    typedef Map<Eigen::MatrixXd> MapMatd;
    
    const int maxiter(as<int>(maxit));
    const double toler(as<double>(tol));
    const Eigen::Map<MatrixXd> XX(as<MapMatd>(X));
    const Eigen::Map<MatrixXd> SS(as<MapMatd>(S));
    const Eigen::Map<VectorXd> YY(as<MapVecd>(Y));
    
    const int n(XX.rows());
    const int p(XX.cols());
    
    MatrixXd SX(MatrixXd(n, p));
    
    SX = SS.llt().solve(XX);
    
    MatrixXd XSX((XX.adjoint()) * SX);
    VectorXd XSY((SX.adjoint()) * YY);
    
    VectorXd beta = conjugate_gradient(XSX, XSY, maxiter, toler);
    
    
    return List::create(Named("beta") = beta);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}


RcppExport SEXP gls_direct(SEXP X, SEXP S, SEXP Y, SEXP maxit, SEXP tol)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    
    using Rcpp::List;
    
    typedef Map<VectorXd> MapVecd;
    typedef Map<Eigen::MatrixXd> MapMatd;
    
    const int maxiter(as<int>(maxit));
    const double toler(as<double>(tol));
    const Eigen::Map<MatrixXd> XX(as<MapMatd>(X));
    const Eigen::Map<MatrixXd> SS(as<MapMatd>(S));
    const Eigen::Map<VectorXd> YY(as<MapVecd>(Y));
    
    const int n(XX.rows());
    const int p(XX.cols());
    
    MatrixXd SX(MatrixXd(n, p));
    
    SX = SS.llt().solve(XX);
    
    MatrixXd XSX((XX.adjoint()) * SX);
    VectorXd XSY((SX.adjoint()) * YY);
    
    VectorXd beta = XSX.ldlt().solve(XSY);
    
    
    return List::create(Named("beta") = beta);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

