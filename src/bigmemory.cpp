#include <iostream>
#include "linop.h"
#include "utils.h"
#include "math.h"
#include <bigmemory/MatrixAccessor.hpp>
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

// Logic for BigColSums.
template <typename T>
NumericVector BigColSums(XPtr<BigMatrix> pMat, MatrixAccessor<T> mat) {
  
  // Create the vector we'll store the column sums in.
  NumericVector colSums(pMat->ncol());
  for (size_t i=0; i < pMat->ncol(); ++i)
    colSums[i] = std::accumulate(mat[i], mat[i]+pMat->nrow(), 0.0);
  return colSums;
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
      Map<MatrixXchar> bM = Map<MatrixXchar>((char *)bMPtr->matrix(), bMPtr->nrow(), bMPtr->ncol()  );
      int p = bM.cols();
      MatrixXchar crossprod = MatrixXchar(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate( bM.adjoint() );
      return wrap(crossprod);
    } else if (type == 2) 
    {
      Map<MatrixXshort> bM = Map<MatrixXshort>((short *)bMPtr->matrix(), bMPtr->nrow(), bMPtr->ncol()  );
      int p = bM.cols();
      MatrixXshort crossprod = MatrixXshort(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate( bM.adjoint() );
      return wrap(crossprod);
    } else if (type == 4) 
    {
      Map<MatrixXi> bM = Map<MatrixXi>((int *)bMPtr->matrix(), bMPtr->nrow(), bMPtr->ncol()  );
      int p = bM.cols();
      MatrixXi crossprod = MatrixXi(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate( bM.adjoint() );
      return wrap(crossprod);
    } else if (type == 8) 
    {
      Map<MatrixXd> bM = Map<MatrixXd>((double *)bMPtr->matrix(), bMPtr->nrow(), bMPtr->ncol()  );
      int p = bM.cols();
      MatrixXd crossprod = MatrixXd(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate( bM.adjoint() );
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
      Map<MatrixXchar> bM = Map<MatrixXchar>((char *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol()  );
      Vectorchar colSums = bM.colwise().sum();
      return wrap(colSums);
    } else if (type == 2) 
    {
      Map<MatrixXshort> bM = Map<MatrixXshort>((short *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol()  );
      Vectorshort colSums = bM.colwise().sum();
      return wrap(colSums);
    } else if (type == 4) 
    {
      Map<MatrixXi> bM = Map<MatrixXi>((int *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol()  );
      VectorXi colSums = bM.colwise().sum();
      return wrap(colSums);
    } else if (type == 8) 
    {
      Map<MatrixXd> bM = Map<MatrixXd>((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol()  );
      VectorXd colSums = bM.colwise().sum();
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



// [[Rcpp::export]]
NumericVector BigColSums(SEXP pBigMat) {
  // First we have to tell Rcpp what class to use for big.matrix objects.
  // This object stores the attributes of the big.matrix object passed to it
  // by R.
  XPtr<BigMatrix> xpMat(pBigMat);
  
  // To access values in the big.matrix, we need to create a MatrixAccessor
  // object of the appropriate type. Note that in every case we are still
  // returning a NumericVector: this is because big.matrix objects only store
  // numeric values in R, even if their type is set to 'char'. The types
  // simply correspond to the number of bytes used for each element.
  switch(xpMat->matrix_type()) {
  case 1:
    return BigColSums(xpMat, MatrixAccessor<char>(*xpMat));
  case 2:
    return BigColSums(xpMat, MatrixAccessor<short>(*xpMat));
  case 4:
    return BigColSums(xpMat, MatrixAccessor<int>(*xpMat));
  case 8:
    return BigColSums(xpMat, MatrixAccessor<double>(*xpMat));
  default:
    // This case should never be encountered unless the implementation of
    // big.matrix changes, but is necessary to implement shut up compiler
    // warnings.
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}




