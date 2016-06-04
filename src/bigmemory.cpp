#include <iostream>
#include "linop.h"
#include "utils.h"
#include "math.h"
#include <bigmemory/BigMatrix.h>

using Eigen::MatrixXf;
using Eigen::VectorXf;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;
using namespace Rcpp;
using namespace RcppEigen;

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> BigEigenColSums(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& bigMat) 
{
  return bigMat.colwise().sum();
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> BigEigenCrossprod(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& bigMat) 
{
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXTT;
  int p = bigMat.cols();
  
  MatrixXTT retmat(p, p);
  retmat.setZero();
  
  return retmat.template selfadjointView<Eigen::Upper>().rankUpdate( bigMat.adjoint() );
}

RcppExport SEXP crossprod_big(SEXP X_)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    
    typedef Eigen::Matrix<char, Eigen::Dynamic, Eigen::Dynamic> MatrixXchar;
    typedef Eigen::Matrix<short, Eigen::Dynamic, Eigen::Dynamic> MatrixXshort;
    
    XPtr<BigMatrix> bMPtr(X_);
    
    
    unsigned int type = bMPtr->matrix_type();
    
    
    if (type == 1) 
    {
      MatrixXchar bM = Map<MatrixXchar>((char *)bMPtr->matrix(), bMPtr->nrow(), bMPtr->ncol()  );
      MatrixXchar crossprod = BigEigenCrossprod(bM);
      return wrap(crossprod);
    } else if (type == 2) 
    {
      MatrixXshort bM = Map<MatrixXshort>((short *)bMPtr->matrix(), bMPtr->nrow(), bMPtr->ncol()  );
      MatrixXshort crossprod = BigEigenCrossprod(bM);
      return wrap(crossprod);
    } else if (type == 4) 
    {
      MatrixXi bM = Map<MatrixXi>((int *)bMPtr->matrix(), bMPtr->nrow(), bMPtr->ncol()  );
      MatrixXi crossprod = BigEigenCrossprod(bM);
      return wrap(crossprod);
    } else if (type == 8) 
    {
      MatrixXd bM = Map<MatrixXd>((double *)bMPtr->matrix(), bMPtr->nrow(), bMPtr->ncol()  );
      MatrixXd crossprod = BigEigenCrossprod(bM);
      return wrap(crossprod);
    } else {
      // We should never get here, but it resolves compiler warnings.
      throw Rcpp::exception("Undefined type for provided big.matrix");
    }
    
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

RcppExport SEXP colsums_big(SEXP X_)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    
    typedef Eigen::Matrix<char, Eigen::Dynamic, Eigen::Dynamic> MatrixXchar;
    typedef Eigen::Matrix<short, Eigen::Dynamic, Eigen::Dynamic> MatrixXshort;
    typedef Eigen::Matrix<char, Eigen::Dynamic, 1> Vectorchar;
    typedef Eigen::Matrix<short, Eigen::Dynamic, 1> Vectorshort;
    
    XPtr<BigMatrix> xpMat(X_);
    
    
    unsigned int type = xpMat->matrix_type();
    // The data stored in the big.matrix can either be represent by 1, 2,
    // 4, or 8 bytes. See the "type" argument in `?big.matrix`.
    if (type == 1) 
    {
      MatrixXchar bM = Map<MatrixXchar>((char *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol()  );
      Vectorchar colSums = BigEigenColSums(bM);
      return wrap(colSums);
    } else if (type == 2) 
    {
      MatrixXshort bM = Map<MatrixXshort>((short *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol()  );
      Vectorshort colSums = BigEigenColSums(bM);
      return wrap(colSums);
    } else if (type == 4) 
    {
      MatrixXi bM = Map<MatrixXi>((int *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol()  );
      VectorXi colSums = BigEigenColSums(bM);
      return wrap(colSums);
    } else if (type == 8) 
    {
      MatrixXd bM = Map<MatrixXd>((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol()  );
      VectorXd colSums = BigEigenColSums(bM);
      return wrap(colSums);
    } else {
      // We should never get here, but it resolves compiler warnings.
      throw Rcpp::exception("Undefined type for provided big.matrix");
    }

  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}
