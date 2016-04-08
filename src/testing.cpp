
#include <igl/slice.h>
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


//port faster cross product 
RcppExport SEXP crossprodsubset(SEXP X, SEXP idxvec_)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::MatrixXi;
    using Eigen::VectorXi;
    using Eigen::Lower;
    typedef float Scalar;
    typedef double Double;
    typedef Eigen::Matrix<Double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Double, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Map<const Matrix> MapMat;
    typedef Eigen::Map<const Vector> MapVec;
    typedef Map<VectorXd> MapVecd;
    typedef Map<VectorXi> MapVeci;
    typedef Eigen::SparseVector<double> SparseVector;
    typedef Eigen::SparseVector<int> SparseVectori;
    
    //const Map<MatrixXd> A(as<Map<MatrixXd> >(X));
    
    const MatrixXd A(as<MatrixXd>(X));
    const VectorXi idxvec(as<VectorXi>(idxvec_));
    //MapMat A(as<MapMat>(X));
    
    int nvars = A.cols() - 1;
    MatrixXd sub;
    
    VectorXi all;
    igl::colon<int>(0, int(nvars), all);
    
    igl::slice(A, idxvec, all, sub);
    
    const int n(A.cols());
    MatrixXd AtA(MatrixXd(n, n).setZero().
                   selfadjointView<Lower>().rankUpdate(sub.adjoint() ));
    
    
    return wrap(AtA);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}


//port faster cross product 
RcppExport SEXP crossprodxval(SEXP X, SEXP idxvec_)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXi;
    using Eigen::Lower;
    using Rcpp::List;
    typedef float Scalar;
    typedef double Double;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixRXd;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Map<Matrix> MapMat;
    typedef Eigen::Map<MatrixRXd> MapRMat;
    typedef Eigen::Map<const Vector> MapVec;
    typedef Map<VectorXd> MapVecd;
    typedef Map<VectorXi> MapVeci;
    typedef Eigen::SparseVector<double> SparseVector;
    typedef Eigen::SparseVector<int> SparseVectori;
    
    
    const MapMat AA(as<MapMat>(X));
    const MapVeci idxvec(as<MapVeci>(idxvec_));
    
    const int n(AA.cols());
    MatrixXd AtA(MatrixXd(n, n));
    MatrixRXd A = AA;
    int nfolds = idxvec.maxCoeff();
    
    List xtxlist(nfolds);
    
    for (int k = 1; k < nfolds + 1; ++k)
    {
    
      VectorXi idxbool = (idxvec.array() == k).cast<int>();
      int nrow = idxbool.size();
      int numelem = idxbool.sum();
      VectorXi idx(numelem);
      int c = 0;
      for (int i = 0; i < nrow; ++i)
      {
        if (idxbool(i) == 1)
        {
          idx(c) = i;
          ++c;
        }
      }
      
      int nvars = A.cols() - 1;
      MatrixXd sub(numelem, n);
      
      for (int r = 0; r < numelem; ++r)
      {
        sub.row(r) = A.row(idx(r));
      }
      
      MatrixXd AtAtmp(MatrixXd(n, n).setZero().
                      selfadjointView<Lower>().rankUpdate(sub.adjoint() ));
      xtxlist[k-1] = AtAtmp;
    }
    
    return wrap(xtxlist);
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





RcppExport SEXP appendMat(SEXP X_, SEXP Y_)
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
    
    const int nn(X.rows());
    const int pp(X.cols());
    
    MatrixXd datX(nn, pp+1);
    
    datX << Y, X ;
    
    return wrap(datX);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}



RcppExport SEXP strcompare(SEXP str_)
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
    
    const std::vector<std::string> strr(as< std::vector<std::string> >(str_));
    
    bool res = strr[0] == "truth";
    int res2 = strr.size();
    
    return wrap(res2);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}


