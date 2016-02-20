#include <iostream>
#include "testing.h"
#include "utils.h"
#include "math.h"

using namespace Rcpp;
using namespace RcppEigen;

//port faster cross product 
RcppExport SEXP crossprodeig(SEXP X)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::Lower;
    typedef float Scalar;
    typedef double Double;
    typedef Eigen::Matrix<Double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Double, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Map<const Matrix> MapMat;
    typedef Eigen::Map<const Vector> MapVec;
    const Eigen::Map<MatrixXd> A(as<Map<MatrixXd> >(X));
    //MapMat A(as<MapMat>(X));
    
    const int n(A.cols());
    MatrixXd AtA(MatrixXd(n, n).setZero().
                   selfadjointView<Lower>().rankUpdate(A.adjoint()));
    
    //std::cout << "crossprod: " << AtA << std::endl;
      
    //Matrix AtAf = Map<Matrix>( xtx_ptr, AtA.rows(), AtA.cols() );
    
    //std::cout << "crossprod2: " << AtAf << std::endl;
    
    MatOpSymLowerDouble<Double> op(AtA);
    Spectra::SymEigsSolver< Double, Spectra::BOTH_ENDS, MatOpSymLowerDouble<Double> > eigs(&op, 2, 5);
    srand(0);
    eigs.init();
    eigs.compute(500, 0.005);
    Vector evals = eigs.eigenvalues();
    
    return wrap(evals);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}



RcppExport SEXP matveccrossprodidx(SEXP X_, SEXP Y_, SEXP idx_)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::VectorXi;
    using Rcpp::List;
    using Eigen::MappedSparseMatrix;
    using Eigen::SparseMatrix;
    using Eigen::Upper;
    typedef MappedSparseMatrix<double> MSpMat;
    typedef SparseMatrix<double> SpMat;
    typedef Map<VectorXi> MapVeci;
    typedef Map<VectorXd> MapVecd;
    typedef Map<MatrixXd> MapMatd;
    
    const MapMatd X(as<MapMatd>(X_));
    const MapVecd Y(as<MapVecd>(Y_));
    const MapVeci idx(as<MapVeci>(idx_));
    //const MapMatd init_(as<MapMatd>(init));
    
    const int nn(X.rows());
    const int pp(X.cols());
    const int rr(idx.size());
    VectorXd retvec(rr);
    
    for (int cl = 0; cl < rr; ++cl)
    {
      //auto jth_column = X.data() + nn * cl;
      //retvec(cl) = (X.col(cl).array() * Y(cl)).sum();
      retvec(cl) = X.col(idx(cl)-1).dot(Y);
      //retvec(cl) = (X.col(idx(cl)-1).array() * Y.array()).sum();
    }
    
    //std::cout << "crossprod: " << AtA << std::endl;
    
    //Matrix AtAf = Map<Matrix>( xtx_ptr, AtA.rows(), AtA.cols() );
    
    
    return wrap(retvec);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}


RcppExport SEXP matvecprodidx(SEXP A_, SEXP b_, SEXP idx_)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::VectorXi;
    using Rcpp::List;
    using Eigen::MappedSparseMatrix;
    using Eigen::SparseMatrix;
    using Eigen::Upper;
    typedef MappedSparseMatrix<double> MSpMat;
    typedef SparseMatrix<double> SpMat;
    typedef Map<VectorXi> MapVeci;
    typedef Map<VectorXd> MapVecd;
    typedef Map<MatrixXd> MapMatd;
    
    const MapMatd A(as<MapMatd>(A_));
    const MapVecd b(as<MapVecd>(b_));
    const MapVeci idx(as<MapVeci>(idx_));
    //const MapMatd init_(as<MapMatd>(init));
    
    const int nn(A.rows());
    const int pp(A.cols());
    const int rr(idx.size());
    VectorXd retvec(nn);
    retvec.setZero();
    
    
    for (int cl = 0; cl < rr; ++cl)
    {
      for (int r = 0; r < nn; ++r)
      {
        retvec(r) += A(r, idx(cl) - 1) * b( idx(cl) - 1 );
      }
    }
    
    /*
    for (size_t cl = 0; cl < rr; ++cl)
      for (size_t r = 0; r < nn; ++r)
      {
        retvec(r) += A(r,idx(cl)-1) * b(idx(cl)-1);
      }
     */
    
    
    return wrap(retvec);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}