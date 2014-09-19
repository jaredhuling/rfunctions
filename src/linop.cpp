
#include "linop.h"
#include "utils.h"

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

inline Eigen::ArrayXd Dplus(const Eigen::ArrayXd& d) {
  Eigen::ArrayXd di(d.size());
  double comp(d.maxCoeff() * std::numeric_limits<double>::epsilon() * d.size());
  for (int j = 0; j < d.size(); ++j) di[j] = (d[j] < comp) ? 0. : 1./d[j];
  return di;
}


//port fast matrix rank function
RcppExport SEXP fastRank(SEXP AA)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::Lower;
    using Eigen::Upper;
    using Eigen::VectorXd;

    typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
    const MapMatd A(as<MapMatd>(AA));
    
    const Eigen::SelfAdjointEigenSolver<MatrixXd> VLV(xtx(A).selfadjointView<Lower>());
    const Eigen::ArrayXd Dp(Dplus(VLV.eigenvalues()).sqrt());
    const int r((Dp > 0).count());
    
    return wrap(r);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}


/////////////////////////////////////////////
//  Solvers for linear systems
/////////////////////////////////////////////

//Biconjugate gradient stabilized method
RcppExport SEXP BiCGSTAB_eigen(SEXP A, SEXP b, SEXP maxit, SEXP tol)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::BiCGSTAB;
    using Rcpp::List;
    typedef Map<VectorXd> MapVecd;
    typedef Map<Eigen::MatrixXd> MapMatd;
    
    const int maxiter(as<int>(maxit));
    const double toler(as<double>(tol));
    Eigen::Map<MatrixXd> AA(as<MapMatd>(A));
    Eigen::Map<VectorXd> bb(as<MapVecd>(b));
    
    BiCGSTAB<MatrixXd > solver;
     
    const int n(AA.cols());
    VectorXd solution(n);
    
    solver.setMaxIterations(maxiter);  
    solver.setTolerance(toler);
    
    solver.compute(AA);
    solution = solver.solve(bb);
    
    return List::create(Named("x") = solution,
                        Named("iters") = solver.iterations(),
                        Named("error") = solver.error());
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}


//Sparse Biconjugate gradient stabilized method
RcppExport SEXP BiCGSTAB_sparse_eigen(SEXP A, SEXP b, SEXP maxit, SEXP tol)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::BiCGSTAB;
    using Eigen::MappedSparseMatrix;
    using Eigen::SparseMatrix;
    using Rcpp::List;
    typedef MappedSparseMatrix<double> MSpMat;
    typedef SparseMatrix<double> SpMat;
    typedef Map<VectorXd> MapVecd;
    typedef Map<Eigen::MatrixXd> MapMatd;
    
    const int maxiter(as<int>(maxit));
    const double toler(as<double>(tol));
    SpMat AA(as<MSpMat>(A));
    MapVecd bb(as<MapVecd>(b));
    
    BiCGSTAB<SpMat > solver;
     
    const int n(AA.cols());
    VectorXd solution(n);
    
    solver.setMaxIterations(maxiter);  
    solver.setTolerance(toler);
    
    solver.compute(AA);
    solution = solver.solve(bb);
    
    return List::create(Named("x") = solution,
                        Named("iters") = solver.iterations(),
                        Named("error") = solver.error());
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}


// Conjugate gradient
RcppExport SEXP conjugate_gradient_eigen(SEXP A, SEXP b, SEXP maxit, SEXP tol)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::ConjugateGradient;
    using Rcpp::List;
    typedef Map<VectorXd> MapVecd;
    typedef Map<Eigen::MatrixXd> MapMatd;
    
    const int maxiter(as<int>(maxit));
    const double toler(as<double>(tol));
    Eigen::Map<MatrixXd> AA(as<MapMatd>(A));
    Eigen::Map<VectorXd> bb(as<MapVecd>(b));
    
    ConjugateGradient<MatrixXd > solver;
     
    const int n(AA.cols());
    VectorXd solution(n);
    
    solver.setMaxIterations(maxiter);  
    solver.setTolerance(toler);
    
    solver.compute(AA);
    solution = solver.solve(bb);
    
    return List::create(Named("x") = solution,
                        Named("iters") = solver.iterations(),
                        Named("error") = solver.error());
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}


// Sparse conjugate gradient
RcppExport SEXP conjugate_gradient_sparse_eigen(SEXP A, SEXP b, SEXP maxit, SEXP tol)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::ConjugateGradient;
    using Eigen::MappedSparseMatrix;
    using Eigen::SparseMatrix;
    using Rcpp::List;
    typedef MappedSparseMatrix<double> MSpMat;
    typedef SparseMatrix<double> SpMat;
    typedef Map<VectorXd> MapVecd;
    typedef Map<Eigen::MatrixXd> MapMatd;
    
    const int maxiter(as<int>(maxit));
    const double toler(as<double>(tol));
    SpMat AA(as<MSpMat>(A));
    MapVecd bb(as<MapVecd>(b));
    
    ConjugateGradient<SpMat > solver;
     
    const int n(AA.cols());
    VectorXd solution(n);
    
    solver.setMaxIterations(maxiter);  
    solver.setTolerance(toler);
    
    solver.compute(AA);
    solution = solver.solve(bb);
    
    return List::create(Named("x") = solution,
                        Named("iters") = solver.iterations(),
                        Named("error") = solver.error());
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

