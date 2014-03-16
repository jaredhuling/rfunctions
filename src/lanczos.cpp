
#include "lanczos.h"

using namespace Rcpp;
using namespace RcppEigen;


//Lanczos Bidiagonalization
RcppExport SEXP GKLBidiag(SEXP A, SEXP v, SEXP k)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Rcpp::List;
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::JacobiSVD;
    typedef Map<VectorXd> MapVecd;
    typedef Map<Eigen::MatrixXd> MapMatd;
    const MapMatd AA(as<MapMatd>(A));
    const MapVecd Vinit(as<MapVecd>(v));
    const int kk(as<int>(k));
    //const bool reorth(as<bool>(reorthog));
    
    VectorXd v(Vinit);
    VectorXd u(AA.rows());
    VectorXd Uprev(v.size());
    VectorXd Vnew(v.size());
    MatrixXd B(MatrixXd::Zero(kk, kk));
    double d(0);
    //VectorXd Unew(v.size());
    
    v.normalize();
    
    VectorXd beta(kk);
    VectorXd alpha(kk);
    
    for (int i = 0; i < kk; i++) {
      Uprev = u;
      u = AA * v;
      if (i > 0) {
        u -= beta(i - 1) * Uprev;
        //add reorthogonalization
      }
      alpha(i) = u.norm();
      u /= alpha(i);
      Vnew = AA.adjoint() * u - alpha(i) * v;
      beta(i) = Vnew.norm();
      v = Vnew.array() / beta(i);
    }
    
    B.diagonal() = alpha;
    for (int i = 0; i < kk - 1; i++) {
      B(i, i + 1) = beta(i);
    }
    
    JacobiSVD<MatrixXd> svd(B);
    d = svd.singularValues()(0);
    
    return List::create(Named("d") = d,
                        Named("u") = u,
                        Named("v") = v,
                        Named("alpha") = alpha,
                        Named("beta") = beta);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

//Lanczos Bidiagonalization for Sparse Matrices
RcppExport SEXP GKLBidiagSparse(SEXP A, SEXP v, SEXP k)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Rcpp::List;
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::MappedSparseMatrix;
    using Eigen::SparseMatrix;
    using Eigen::VectorXd;
    using Eigen::JacobiSVD;
    typedef MappedSparseMatrix<double> MSpMat;
    typedef SparseMatrix<double> SpMat;
    typedef Map<VectorXd> MapVecd;
    typedef Map<MatrixXd> MapMatd;
    
    const SpMat AA(as<MSpMat>(A));
    const MapVecd Vinit(as<MapVecd>(v));
    const int kk(as<int>(k));
    //const bool reorth(as<bool>(reorthog));
    
    VectorXd v(Vinit);
    VectorXd u(AA.rows());
    VectorXd Uprev(v.size());
    VectorXd Vnew(v.size());
    MatrixXd B(MatrixXd::Zero(kk, kk));
    double d(0);
    //VectorXd Unew(v.size());
    
    v.normalize();
        
    VectorXd beta(kk);
    VectorXd alpha(kk);
    
    for (int i = 0; i < kk; i++) {
      Uprev = u;
      u = AA * v;
      if (i > 0) {
        u -= beta(i - 1) * Uprev;
        //add reorthogonalization
      }
      alpha(i) = u.norm();
      u /= alpha(i);
      Vnew = AA.adjoint() * u - alpha(i) * v;
      beta(i) = Vnew.norm();
      v = Vnew.array() / beta(i);
    }
    
    B.diagonal() = alpha;
    for (int i = 0; i < kk - 1; i++) {
      B(i, i + 1) = beta(i);
    }
    
    JacobiSVD<MatrixXd> svd(B);
    d = svd.singularValues()(0);
    
    return List::create(Named("d") = d,
                        Named("u") = u,
                        Named("v") = v,
                        Named("alpha") = alpha,
                        Named("beta") = beta);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

//port faster cross product 
RcppExport SEXP BidiagPoly(SEXP X, SEXP alpha, SEXP beta)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    typedef Eigen::Map<VectorXd> MapVecd;
    const MapVecd alph(as<MapVecd>(alpha));
    const MapVecd bet(as<MapVecd>(beta));
    const double xx(as<double>(X));
    
    double p0(1 / alph(0));
    double q0(1);
    double p(p0);
    double q(q0);
    
    for (int k = 0; k < alph.size() - 1; k++) {
      q = (xx * p0 - alph(k) * q0) / bet(k);
      p  = (q - bet(k) * p0) / alph(k + 1);
      p0 = p;
      q0 = q;
    }
    
    return wrap(p);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}
