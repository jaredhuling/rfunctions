

#include "multivar_norm.h"


RcppExport SEXP multivarNorm(SEXP n, SEXP mean, SEXP sigma, SEXP seed) 
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::Map;
    using Rcpp::as;
    typedef Eigen::Map<VectorXd> MapVecd;
    typedef Map<Eigen::MatrixXd> MapMatd;

    const MapVecd mean_(as<MapVecd>(mean));
    const Eigen::Map<MatrixXd> sigma_(as<MapMatd>(sigma));
    const int n_(as<int>(n));
    const int seed_(as<int>(seed));
    
    int size = sigma_.cols();
        
    Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
    
    if (seed_ != -999999) {
      Eigen::internal::scalar_normal_dist_op<double>::rng.seed(seed_); // Seed the rng
    }
    
    MatrixXd normTransform(size,size);

    Eigen::LLT<MatrixXd> cholSolver(sigma_);

    // We can only use the cholesky decomposition if 
    // the covariance matrix is symmetric, pos-definite.
    // But a covariance matrix might be pos-semi-definite.
    // In that case, we'll go to an EigenSolver
    if (cholSolver.info()==Eigen::Success) {
      // Use cholesky solver
      normTransform = cholSolver.matrixL();
    } else {
      // Use eigen solver
      Eigen::SelfAdjointEigenSolver<MatrixXd> eigenSolver(sigma_);
      normTransform = eigenSolver.eigenvectors() 
                      * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }

    MatrixXd samples = (normTransform 
                        * MatrixXd::NullaryExpr(size, n_, randN)).colwise() 
                        + mean_;

    return wrap(samples.adjoint());
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}