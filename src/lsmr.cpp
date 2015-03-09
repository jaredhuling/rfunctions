#include <iostream>
#include "lsmr.h"
#include "utils.h"
#include "math.h"

using namespace Rcpp;
using namespace RcppEigen;

//port faster cross product 
RcppExport SEXP lsmr_sparse(SEXP A,        SEXP b,
                            SEXP lambda, 
                            SEXP atol,     SEXP btol, 
                            SEXP conlim,   SEXP maxit, 
                            SEXP localSize)
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
    
    
    int maxiter(as<int>(maxit));
    const int locSize(as<int>(localSize));
    const double lam(as<double>(lambda));
    const double btoler(as<double>(btol));
    const double atoler(as<double>(atol));
    double clim(as<double>(conlim));
    const SpMat AA(as<MSpMat>(A));
    const MapMatd bb(as<MapMatd>(b));
    
    const int n(AA.cols());
    const int m(AA.rows());
    
    VectorXd u(m);
    VectorXd v(n);
    MatrixXd localV;
    
    // Determine dimensions m and n, and
    // form the first vectors u and v.
    // These satisfy  beta*u = b,  alpha*v = A'u.
    
    u = bb;
    double beta(u.norm());
    if (beta > 0) {
      u = u.array() / beta;
    }
    
    v = AA.adjoint() * u;
    
    int mindim = std::min(m, n);
    
    maxiter = std::min(maxiter, mindim);
    
    double alpha(v.norm());
    if (alpha > 0) {
      v = v.array() / alpha;
    }
    
    // Initialization for local reorthogonalization.
    int localPointer;
    bool localOrtho = false;
    bool localVQueueFull = true;
    if (locSize > 0) {
      localPointer = 0;
      localOrtho = true;
      localVQueueFull = false;
      
      // Preallocate storage for the relevant number of latest v_k's.
      localV = MatrixXd::Zero(n, std::min(locSize, mindim));
    }
    
    // Initialize variables for 1st iteration.
    double zetabar = alpha * beta;
    double alphabar = alpha;
    double rho = 1;
    double rhobar = 1;
    double cbar = 1;
    double sbar = 0;
    double condA = 1;
    int istop;
    
    VectorXd h(v);
    
    VectorXd hbar(MatrixXd::Zero(n, 1));
    VectorXd x(MatrixXd::Zero(n, 1));
    
    // Initialize variables for estimation of ||r||.
    double betadd = beta;
    double betad = 0;
    double betaacute = 0;
    double betacheck = 0;
    double rhodold = 1;
    double tautildeold = 0;
    double thetatilde = 0;
    double zeta = 0;
    double d = 0;
    
    double thetatildeold = 0;
    double rhotildeold   = 0;
    double ctildeold     = 0;
    double stildeold     = 0;
    
    double alphahat = 0;
    double chat = 0;
    double shat = 0;
    double rhoold = 0;
    double c = 0;
    double s = 0;
    double normA = 0;
    double thetanew = 0;
    VectorXd vec1(2);
    
    double rhobarold = 0;
    double zetaold = 0;
    double rhotemp = 0;
    double thetabar;
    double taud;
    double betahat;
    
    
    // Initialize variables for estimation of ||A|| and cond(A).
    double normA2 = pow(alpha, 2);
    double normx = 1;
    double maxrbar = 0;
    double minrbar = 1e+100;
    
    // Items for use in stopping rules.
    double normb  = beta;
    double test1 = 1;
    double test2 = 1;
    double test3 = 1;
    double t1 = 1;
    double rtoler = 1;
    double ctol   = 0;         
    if (clim > 0) {
      ctol = 1 / clim;
    }
    double normr  = beta;
    
    // Exit if b=0 or A'b = 0.
    double normAr = alpha * beta;
    
    int itn;
    for (itn = 0; itn < maxiter; itn++) {
      
      // Perform the next step of the bidiagonalization to obtain the
      // next beta, u, alpha, v.  These satisfy the relations
      // beta*u  =  A*v  - alpha*u,
      // alpha*v  =  A'*u - beta*v.
      
      u = (AA * v).array() - (alpha * u.array()).array();
      
      beta = u.norm();
      normb = beta;
      
      if (beta > 0) {
        u = u.array() / beta;
        //if (localOrtho) {
        //  localVEnqueue(v);
        //}
        
        v = (AA.adjoint() * u).array() - (beta * v.array()).array();
        
        //if (localOrtho) {
        //  v = localVOrtho(v);
        //}
        
        alpha = v.norm();
        if (alpha > 0) {  
          v = v.array() / alpha;
        }       
        
      }
      
      // At this point, beta = beta_{k+1}, alpha = alpha_{k+1}.
  
      // Construct rotation Qhat_{k,2k+1}.
      vec1 << alphabar, lam;
      alphahat = vec1.norm();
      chat     = alphabar / alphahat;
      shat     = lam / alphahat;
        
        
      // Use a plane rotation (Q_i) to turn B_i to R_i.
      rhoold   = rho;
      vec1 << alphahat, beta;
      rho      = vec1.norm();
      c        = alphahat / rho;
      s        = beta / rho;
      thetanew = s * alpha;
      alphabar = c * alpha;
  
      // Use a plane rotation (Qbar_i) to turn R_i^T to R_i^bar.
      rhobarold = rhobar;
      zetaold   = zeta;
      thetabar  = sbar * rho;
      rhotemp   = cbar * rho;
      vec1 << cbar*rho, thetanew;
      rhobar    = vec1.norm();
      cbar      = (cbar * rho) / rhobar;
      sbar      = thetanew / rhobar;
      zeta      =   cbar * zetabar;
      zetabar   = - sbar * zetabar;
  
      // Update h, h_hat, x.
      hbar      = h.array() - (thetabar * rho / (rhoold * rhobarold)) * hbar.array();
      x         = x.array() + (zeta / (rho*rhobar)) * hbar.array();
      h         = v.array() - (thetanew / rho) * h.array();
      
      // Estimate of ||r||.
    
      // Apply rotation Qhat_{k,2k+1}.
      betaacute =   chat * betadd;
      betacheck = - shat * betadd;
  
      // Apply rotation Q_{k,k+1}.
      betahat   =   c * betaacute;
      betadd    = - s * betaacute;
        
      // Apply rotation Qtilde_{k-1}.
      // betad = betad_{k-1} here.
  
      thetatildeold = thetatilde;
      vec1 << rhodold, thetabar;
      rhotildeold   = vec1.norm();
      ctildeold     = rhodold / rhotildeold;
      stildeold     = thetabar / rhotildeold;
      thetatilde    = stildeold * rhobar;
      rhodold       =   ctildeold * rhobar;
      betad         = - stildeold * betad + ctildeold * betahat;
  
      // betad   = betad_k here.
      // rhodold = rhod_k  here.
  
      tautildeold   = (zetaold - thetatildeold * tautildeold) / rhotildeold;
      taud          = (zeta - thetatilde * tautildeold) / rhodold;
      d             = d + pow(betacheck, 2);
      normr         = sqrt(d + pow((betad - taud), 2) + pow(betadd, 2) );
      
      // Estimate ||A||.
      normA2        = normA2 + pow(beta, 2);
      normA         = sqrt(normA2);
      normA2        = normA2 + pow(alpha, 2);
      
      // Estimate cond(A).
      maxrbar       = std::max(maxrbar,rhobarold);
      if (itn > 1) { 
        minrbar     = std::min(minrbar,rhobarold);
      }
      condA         = std::max(maxrbar,rhotemp) / std::min(minrbar,rhotemp);

      // Test for convergence.

      // Compute norms for convergence testing.
      normAr  = std::abs(zetabar);
      normx   = x.norm();

      // Now use these norms to estimate certain other quantities,
      // some of which will be small near a solution.
      
      test1   = normr / normb;
      test2   = normAr / (normA * normr);
      test3   = 1 / condA;
      t1      = test1 / (1 + normA * normx / normb);
      rtoler    = btoler + atoler * normA * normx / normb;
      
      // The following tests guard against extremely small values of
      // atol, btol or ctol.  (The user may have set any or all of
      // the parameters atol, btol, conlim  to 0.)
      // The effect is equivalent to the normAl tests using
      // atol = eps,  btol = eps,  conlim = 1/eps.
  
      if (itn >= maxiter) {
        istop = 7; 
      }
      if (1 + test3  <= 1) { 
        istop = 6;
      }
      if (1 + test2  <= 1) {
        istop = 5;
      }
      if (1 + t1     <= 1) {
        istop = 4; 
      }
  
      // Allow for tolerances set by the user.
  
      if  (test3 <= ctol) {
        istop = 3; 
      }
      if  (test2 <= atoler) {
        istop = 2; 
      }
      if  (test1 <= rtoler)  {
        istop = 1; 
      }
      
      if (istop > 0) {
        break;
      }

      
    }
    


    return List::create(Named("x") = x,
                        Named("istop") = istop,
                        Named("iters") = itn,
                        Named("normr") = normr,
                        Named("normAr") = normAr,
                        Named("normA") = normA,
                        Named("condA") = condA,
                        Named("normx") = normx);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}