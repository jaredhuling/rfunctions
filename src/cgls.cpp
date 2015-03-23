#include <iostream>
#include "cgls.h"
#include "utils.h"
#include "math.h"

using namespace Rcpp;
using namespace RcppEigen;


// conjugate gradient least squares
RcppExport SEXP cgls(SEXP A, SEXP b, SEXP lambda, SEXP maxit, SEXP tol, SEXP init)
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
    const double lam(as<double>(lambda));
    const Eigen::Map<MatrixXd> AA(as<MapMatd>(A));
    const Eigen::Map<VectorXd> bb(as<MapVecd>(b));
    const Eigen::Map<VectorXd> xinit(as<MapVecd>(init));
    
    const int n(AA.cols());
    const int m(AA.rows());
    VectorXd x(n);
    VectorXd r(m);
    VectorXd p(n);
    VectorXd s(n);
    VectorXd q(m);
    x = xinit;
    
    double norms0;
    double norms;
    double gamma;
    double gammaold;
    double delta;
    double normx;
    double alpha;
    int iters = maxiter; 
    
    r = bb - AA * x; 
    s = AA.transpose() * r - (lam * x.array()).matrix();
    p = s;
    
    gamma = s.squaredNorm();
    norms0 = sqrt(gamma);
    normx = x.norm();
    
    for (int i = 0; i < maxiter; i++) {
      q = AA * p;
      delta = q.squaredNorm()  +  lam * p.squaredNorm();
      alpha = gamma / delta;
      x = x + (alpha * p.array()).matrix();
      r = r - (alpha * q.array()).matrix();
      
      s = AA.transpose() * r - (lam * x.array()).matrix();
      
      norms = s.norm();
      gammaold = gamma;
      gamma = pow(norms, 2);
      
      normx = x.norm();
      
      if (norms < norms0 * toler || normx * toler >= 1) {
        iters = i;
        break;
      }
      p = s + ((gamma / gammaold) * p.array()).matrix();
    }
    
    
    return List::create(Named("x") = x,
                        Named("iters") = iters);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}


// block conjugate gradient least squares
RcppExport SEXP block_cgls(SEXP A, SEXP b, SEXP lambda, SEXP maxit, SEXP tol, SEXP init)
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
    const double lam(as<double>(lambda));
    const Eigen::Map<MatrixXd> AA(as<MapMatd>(A));
    const Eigen::Map<MatrixXd> bb(as<MapMatd>(b));
    const Eigen::Map<MatrixXd> xinit(as<MapMatd>(init));
    
    const int n(AA.cols());
    const int m(AA.rows());
    const int pp(bb.cols());
    MatrixXd x(MatrixXd(n, pp));
    MatrixXd r(MatrixXd(m, pp));
    MatrixXd p(MatrixXd(n, pp));
    MatrixXd s(MatrixXd(n, pp));
    MatrixXd q(MatrixXd(m, pp));
    x = xinit;
    
    double norms0;
    double norms;
    MatrixXd gamma(MatrixXd(pp, pp));
    MatrixXd gammaold(MatrixXd(pp, pp));
    MatrixXd delta(MatrixXd(pp, pp));
    double normx;
    MatrixXd alpha(MatrixXd(pp, pp));
    int iters = maxiter; 
    
    r = bb - AA * x; 
    s = AA.transpose() * r - (lam * x.array()).matrix();
    p = s;
    
    gamma = xtx(s);
    norms0 = sqrt(gamma.diagonal().sum());
    normx = x.norm();
    
    for (int i = 0; i < maxiter; i++) {
      q = AA * p;
      delta = xtx(q) + (lam * xtx(p).array()).matrix();
      alpha = delta.colPivHouseholderQr().solve(gamma);
      x = x + p * alpha;
      r = r - q * alpha;
      
      s = AA.transpose() * r - (lam * x.array()).matrix();

      gammaold = gamma;
      norms = sqrt(gamma.diagonal().sum());
      gamma = xtx(s);
      
      normx = x.norm();
      
      if (norms < norms0 * toler || normx * toler >= 1) {
        iters = i;
        break;
      }
      p = s + p * gammaold.colPivHouseholderQr().solve(gamma);
    }
    
    
    return List::create(Named("x") = x,
                        Named("iters") = iters);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}



// conjugate gradient least squares
RcppExport SEXP cgls_sparse(SEXP A, SEXP b, SEXP lambda, SEXP maxit, SEXP tol, SEXP init)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Rcpp::List;
    using Eigen::MappedSparseMatrix;
    using Eigen::SparseMatrix;
    
    typedef MappedSparseMatrix<double> MSpMat;
    typedef SparseMatrix<double> SpMat;
    typedef Map<VectorXd> MapVecd;
    typedef Map<Eigen::MatrixXd> MapMatd;
    
    const int maxiter(as<int>(maxit));
    const double toler(as<double>(tol));
    const double lam(as<double>(lambda));
    const SpMat AA(as<MSpMat>(A));
    const Eigen::Map<VectorXd> bb(as<MapVecd>(b));
    const Eigen::Map<VectorXd> xinit(as<MapVecd>(init));
    
    const int n(AA.cols());
    const int m(AA.rows());
    VectorXd x(n);
    VectorXd r(m);
    VectorXd p(n);
    VectorXd s(n);
    VectorXd q(m);
    x = xinit;
    
    double norms0;
    double norms;
    double gamma;
    double gammaold;
    double delta;
    double normx;
    double alpha;
    int iters = maxiter; 
    
    r = bb - AA * x; 
    s = AA.transpose() * r - (lam * x.array()).matrix();
    p = s;
    
    gamma = s.squaredNorm();
    norms0 = sqrt(gamma);
    normx = x.norm();
    
    for (int i = 0; i < maxiter; i++) {
      q = AA * p;
      delta = q.squaredNorm()  +  lam * p.squaredNorm();
      alpha = gamma / delta;
      x = x + (alpha * p.array()).matrix();
      r = r - (alpha * q.array()).matrix();
      
      s = AA.transpose() * r - (lam * x.array()).matrix();
      
      norms = s.norm();
      gammaold = gamma;
      gamma = pow(norms, 2);
      
      normx = x.norm();
      
      if (norms < norms0 * toler || normx * toler >= 1) {
        iters = i;
        break;
      }
      p = s + ((gamma / gammaold) * p.array()).matrix();
    }
    
    
    return List::create(Named("x") = x,
                        Named("iters") = iters);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}


// block conjugate gradient least squares
RcppExport SEXP block_cgls_sparse(SEXP A, SEXP b, SEXP lambda, SEXP maxit, SEXP tol, SEXP init)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Rcpp::List;
    using Eigen::MappedSparseMatrix;
    using Eigen::SparseMatrix;
    
    typedef MappedSparseMatrix<double> MSpMat;
    typedef SparseMatrix<double> SpMat;
    typedef Map<VectorXd> MapVecd;
    typedef Map<Eigen::MatrixXd> MapMatd;
    
    const int maxiter(as<int>(maxit));
    const double toler(as<double>(tol));
    const double lam(as<double>(lambda));
    const SpMat AA(as<MSpMat>(A));
    const Eigen::Map<MatrixXd> bb(as<MapMatd>(b));
    const Eigen::Map<MatrixXd> xinit(as<MapMatd>(init));
    
    const int n(AA.cols());
    const int m(AA.rows());
    const int pp(bb.cols());
    MatrixXd x(MatrixXd(n, pp));
    MatrixXd r(MatrixXd(m, pp));
    MatrixXd p(MatrixXd(n, pp));
    MatrixXd s(MatrixXd(n, pp));
    MatrixXd q(MatrixXd(m, pp));
    x = xinit;
    
    double norms0;
    double norms;
    MatrixXd gamma(MatrixXd(pp, pp));
    MatrixXd gammaold(MatrixXd(pp, pp));
    MatrixXd delta(MatrixXd(pp, pp));
    double normx;
    MatrixXd alpha(MatrixXd(pp, pp));
    int iters = maxiter; 
    
    r = bb - AA * x; 
    s = AA.transpose() * r - (lam * x.array()).matrix();
    p = s;
    
    gamma = xtx(s);
    norms0 = sqrt(gamma.diagonal().sum());
    normx = x.norm();
    
    for (int i = 0; i < maxiter; i++) {
      q = AA * p;
      delta = xtx(q) + (lam * xtx(p).array()).matrix();
      alpha = delta.colPivHouseholderQr().solve(gamma);
      x = x + p * alpha;
      r = r - q * alpha;
      
      s = AA.transpose() * r - (lam * x.array()).matrix();

      gammaold = gamma;
      norms = sqrt(gamma.diagonal().sum());
      gamma = xtx(s);
      
      normx = x.norm();
      
      if (norms < norms0 * toler || normx * toler >= 1) {
        iters = i;
        break;
      }
      p = s + p * gammaold.colPivHouseholderQr().solve(gamma);
    }
    
    
    return List::create(Named("x") = x,
                        Named("iters") = iters);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

