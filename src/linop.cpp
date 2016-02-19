#include <iostream>
#include "linop.h"
#include "utils.h"
#include "math.h"

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
RcppExport SEXP crossprodSparsecpp(SEXP X)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::MappedSparseMatrix;
    using Eigen::SparseMatrix;
    using Eigen::Upper;
    using Eigen::Lower;
    typedef MappedSparseMatrix<double> MSpMat;
    typedef SparseMatrix<double> SpMat;
    const SpMat A(as<MSpMat>(X));
    
    const int n(A.cols());
    SpMat AtA(SpMat(n, n).selfadjointView<Upper>().rankUpdate(A.adjoint()));
    return wrap(AtA);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

//port faster cross product 
RcppExport SEXP crossprodSparsecpphaha(SEXP X, SEXP Y)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::MappedSparseMatrix;
    using Eigen::SparseMatrix;
    using Eigen::Upper;
    using Eigen::Lower;
    typedef MappedSparseMatrix<double> MSpMat;
    typedef SparseMatrix<double> SpMat;
    const SpMat A(as<MSpMat>(X));
    const SpMat B(as<MSpMat>(Y));
    
    const int n(A.cols());
    const int m(B.cols());
    SpMat AtA(n, n);
    SpMat res(n, m);
    AtA.selfadjointView<Lower>().rankUpdate(A.adjoint());
    res = AtA * B;
    return wrap(res);
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
    //using Eigen::IncompleteLUT;
    using Rcpp::List;
    typedef Map<VectorXd> MapVecd;
    typedef Map<Eigen::MatrixXd> MapMatd;
    
    const int maxiter(as<int>(maxit));
    const double toler(as<double>(tol));
    Eigen::Map<MatrixXd> AA(as<MapMatd>(A));
    Eigen::Map<VectorXd> bb(as<MapVecd>(b));
    //BiCGSTAB< MatrixXd,IncompleteLUT<double> > solver;
    BiCGSTAB< MatrixXd > solver;
     
    const int n(AA.cols());
    VectorXd solution(n);
    
    //solver.preconditioner().setFillfactor(7);
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


// Conjugate gradient
RcppExport SEXP conjugate_gradient(SEXP A, SEXP b, SEXP maxit, SEXP tol)
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
    const Eigen::Map<MatrixXd> AA(as<MapMatd>(A));
    const Eigen::Map<VectorXd> bb(as<MapVecd>(b));
    
    const int n(AA.cols());
    VectorXd x(n);
    VectorXd r(n);
    VectorXd p(n);
    VectorXd Ap(n);
    x.fill(0);
    
    double rsold;
    double rsnew;
    double alpha;
    int iters = maxiter; 
    
    r = bb; 
    p = r;
    rsold = r.squaredNorm();
    
    for (int i = 0; i < maxiter; i++) {
      Ap = AA * p;
      alpha = rsold / (p.transpose() * Ap);
      x = x + (alpha * p.array()).matrix();
      r = r - (alpha * Ap.array()).matrix();
      rsnew = r.squaredNorm();
      if (sqrt(rsnew) < toler) {
        iters = i;
        break;
      }
      p = r + ((rsnew / rsold) * p.array()).matrix();
      rsold = rsnew;
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


// Conjugate gradient
RcppExport SEXP ILU_prec_conjugate_gradient_sparse(SEXP A, SEXP b, SEXP maxit, SEXP tol)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
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
    const Eigen::Map<VectorXd> bb(as<MapVecd>(b));
    
    const int n(AA.cols());
    VectorXd x(n);
    VectorXd r(n);
    VectorXd z(n);
    VectorXd p(n);
    VectorXd Ap(n);
    x.fill(0);
    
    Eigen::IncompleteLUT<double> ILUTC;
    ILUTC.setFillfactor(5);
    ILUTC.compute(AA);
    
    double rsold;
    double rsnew;
    double alpha;
    int iters = maxiter; 
    
    r = bb;
    z = ILUTC.solve(r);
    p = r;
    rsold = r.dot(z);

    for (int i = 0; i < maxiter; i++) {
      Ap = AA * p;
      alpha = rsold / (p.transpose() * Ap);
      x = x + (alpha * p.array()).matrix();
      r = r - (alpha * Ap.array()).matrix();
      rsnew = r.squaredNorm();
      if (sqrt(rsnew) < toler) {
        iters = i;
        break;
      }
      z = ILUTC.solve(r);
      rsnew = r.dot(z);
      p = r + ((rsnew / rsold) * p.array()).matrix();
      rsold = rsnew;
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


// Conjugate gradient sparse
RcppExport SEXP conjugate_gradient_sparse(SEXP A, SEXP b, SEXP maxit, SEXP tol)
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
    using Eigen::Upper;
    typedef MappedSparseMatrix<double> MSpMat;
    typedef SparseMatrix<double> SpMat;
    typedef Map<VectorXd> MapVecd;
    typedef Map<MatrixXd> MapMatd;
    
    const int maxiter(as<int>(maxit));
    const double toler(as<double>(tol));
    const SpMat AA(as<MSpMat>(A));
    const Eigen::Map<VectorXd> bb(as<MapVecd>(b));
    
    const int n(AA.cols());
    VectorXd x(n);
    VectorXd r(n);
    VectorXd p(n);
    VectorXd Ap(n);
    x.fill(0);
    
    double rsold;
    double rsnew;
    double alpha;
    int iters = maxiter; 
    
    r = bb; 
    p = r;
    rsold = r.squaredNorm();
    
    for (int i = 0; i < maxiter; i++) {
      Ap = AA * p;
      alpha = rsold / (p.transpose() * Ap);
      x = x + (alpha * p.array()).matrix();
      r = r - (alpha * Ap.array()).matrix();
      rsnew = r.squaredNorm();
      if (sqrt(rsnew) < toler) {
        iters = i;
        break;
      }
      p = r + ((rsnew / rsold) * p.array()).matrix();
      rsold = rsnew;
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


// block Conjugate gradient sparse
RcppExport SEXP block_conjugate_gradient(SEXP A, SEXP b, SEXP maxit, SEXP tol)
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
    using Eigen::Upper;
    typedef MappedSparseMatrix<double> MSpMat;
    typedef SparseMatrix<double> SpMat;
    typedef Map<VectorXd> MapVecd;
    typedef Map<MatrixXd> MapMatd;
    
    const int maxiter(as<int>(maxit));
    const double toler(as<double>(tol));
    const Map<MatrixXd> AA(as<MapMatd>(A));
    const MapMatd bb(as<MapMatd>(b));
    //const MapMatd init_(as<MapMatd>(init));
    
    const int n(AA.cols());
    const int pp(bb.cols());
    MatrixXd x(MatrixXd(n, pp));
    MatrixXd r(MatrixXd(n, pp));
    //MatrixXd theQ(MatrixXd(n, n));
    MatrixXd p(MatrixXd(n, pp));
    MatrixXd Ap(MatrixXd(n, pp));
    x.fill(0);
    
    MatrixXd rsold(MatrixXd(pp, pp));
    MatrixXd rsnew(MatrixXd(pp, pp));
    MatrixXd alpha(MatrixXd(pp, pp));
    int iters = maxiter; 
    
    //MatrixXd RR;


    //Eigen::ColPivHouseholderQR<MatrixXd> pqr(-bb);
    

    
    //int m_r = pqr.rank();

    // Q
    //theQ = pqr.householderQ(); //.setLength(pqr.nonzeroPivots());
    //r = theQ.leftCols(m_r);
    r = -bb;

    
    // R
    //RR = pqr.matrixQR().topRightCorner(m_r, m_r).triangularView<Upper>();

    
    p = -r;
    //rsold = r.colwise().squaredNorm();
    
    rsold = xtx(r);
    

    
    for (int i = 0; i < maxiter; i++) {
      Ap = AA * p;
      alpha = (p.transpose() * Ap).colPivHouseholderQr().solve(rsold);
      x = x + p * alpha;
      r = r + Ap * alpha;
      //rsnew = r.colwise().squaredNorm();
      rsnew = xtx(r);

      if (sqrt(rsnew.diagonal().sum()) < toler) {
        iters = i;
        break;
      }
      p = -r + p * rsold.colPivHouseholderQr().solve(rsnew);
      rsold = rsnew;
    }
    
    
    
    return List::create(Named("x") = x, // x * RR, 
                        Named("iters") = iters);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}


// block Conjugate gradient sparse
RcppExport SEXP block_conjugate_gradient_sparse(SEXP A, SEXP b, SEXP maxit, SEXP tol)
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
    using Eigen::Upper;
    typedef MappedSparseMatrix<double> MSpMat;
    typedef SparseMatrix<double> SpMat;
    typedef Map<VectorXd> MapVecd;
    typedef Map<MatrixXd> MapMatd;
    
    const int maxiter(as<int>(maxit));
    const double toler(as<double>(tol));
    const SpMat AA(as<MSpMat>(A));
    const MapMatd bb(as<MapMatd>(b));
    //const MapMatd init_(as<MapMatd>(init));
    
    const int n(AA.cols());
    const int pp(bb.cols());
    MatrixXd x(MatrixXd(n, pp));
    MatrixXd r(MatrixXd(n, pp));
    //MatrixXd theQ(MatrixXd(n, n));
    MatrixXd p(MatrixXd(n, pp));
    MatrixXd Ap(MatrixXd(n, pp));
    x.fill(0);
    
    MatrixXd rsold(MatrixXd(pp, pp));
    MatrixXd rsnew(MatrixXd(pp, pp));
    MatrixXd alpha(MatrixXd(pp, pp));
    int iters = maxiter; 
    
    //MatrixXd RR;


    //Eigen::ColPivHouseholderQR<MatrixXd> pqr(-bb);
    

    
    //int m_r = pqr.rank();

    // Q
    //theQ = pqr.householderQ(); //.setLength(pqr.nonzeroPivots());
    //r = theQ.leftCols(m_r);
    r = -bb;

    
    // R
    //RR = pqr.matrixQR().topRightCorner(m_r, m_r).triangularView<Upper>();

    
    p = -r;
    //rsold = r.colwise().squaredNorm();
    
    rsold = xtx(r);
    

    
    for (int i = 0; i < maxiter; i++) {
      Ap = AA * p;
      alpha = (p.transpose() * Ap).colPivHouseholderQr().solve(rsold);
      x = x + p * alpha;
      r = r + Ap * alpha;
      //rsnew = r.colwise().squaredNorm();
      rsnew = xtx(r);

      if (sqrt(rsnew.diagonal().sum()) < toler) {
        iters = i;
        break;
      }
      p = -r + p * rsold.colPivHouseholderQr().solve(rsnew);
      rsold = rsnew;
    }
    
    
    
    return List::create(Named("x") = x, // x * RR, 
                        Named("iters") = iters);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}



//port incomplete LU factorization
RcppExport SEXP incompleteLUT_eig(SEXP A, SEXP fillfactor)
{
  using namespace Rcpp;
  using namespace RcppEigen;
  try {
    using Eigen::Map;
    using Eigen::Lower;
    using Eigen::MappedSparseMatrix;
    using Eigen::MatrixXd;
    using Eigen::SparseMatrix;
    typedef Eigen::MappedSparseMatrix<double> MSpMat;
    typedef Eigen::SparseMatrix<double> SpMat;
    const SpMat A_(as<MSpMat>(A));

    
    const int fillfact(as<int>(fillfactor));
    
    Eigen::IncompleteLUT<double> ILUTC;
    ILUTC.setFillfactor(fillfact);
    ILUTC.compute(A_);
    const int numrows(A_.rows());
    const int numcols(A_.cols());
    
    MatrixXd Id(MatrixXd::Identity(numrows, numcols));
    
    return wrap(ILUTC.solve(Id));
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}







RcppExport SEXP sliced_crossprod(SEXP X_, SEXP Y_, SEXP idx_)
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
      retvec(cl) = X.col(idx(cl)-1).dot(Y);
    }
    
    return wrap(retvec);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}







