
#include "linop.h"

using namespace Rcpp;
using namespace RcppEigen;

//port faster cross product 
RcppExport SEXP crossprodcpp(SEXP X)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::Lower;
    const Eigen::Map<MatrixXd> A(as<Map<MatrixXd> >(X));
    const int n(A.cols());
    MatrixXd AtA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(A.adjoint()));
    return wrap(AtA);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

//port faster cross product 
RcppExport SEXP xpwx(SEXP X, SEXP W)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::Lower;
    const Eigen::Map<MatrixXd> A(as<Map<MatrixXd> >(X));
    const Eigen::Map<VectorXd> diag(as<Map<VectorXd> >(W));
    const int n(A.cols());
    MatrixXd AtA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(A.adjoint() * diag.array().sqrt().matrix().asDiagonal()));
    return wrap(AtA);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

//port faster subtract
RcppExport SEXP subcpp(SEXP BB, SEXP CC)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
    const MapMatd B(as<MapMatd>(BB));
    const MapMatd C(as<MapMatd>(CC));
    return wrap(B - C);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}


//port faster add
RcppExport SEXP addcpp(SEXP BB, SEXP CC)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
    const MapMatd B(as<MapMatd>(BB));
    const MapMatd C(as<MapMatd>(CC));
    return wrap(B + C);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}



//port faster subtract sparse matrices
RcppExport SEXP subSparsecpp(SEXP BB, SEXP CC)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::MappedSparseMatrix;
    using Eigen::SparseMatrix;
    typedef Eigen::MappedSparseMatrix<double> MSpMat;
    typedef Eigen::SparseMatrix<double> SpMat;
    const SpMat B(as<MSpMat>(BB));
    const SpMat C(as<MSpMat>(CC));
    return wrap(B - C);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}


//port faster subaddtract sparse matrices
RcppExport SEXP addSparsecpp(SEXP BB, SEXP CC)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::MappedSparseMatrix;
    using Eigen::SparseMatrix;
    typedef Eigen::MappedSparseMatrix<double> MSpMat;
    typedef Eigen::SparseMatrix<double> SpMat;
    const SpMat B(as<MSpMat>(BB));
    const SpMat C(as<MSpMat>(CC));
    return wrap(B + C);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

