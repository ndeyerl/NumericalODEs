/* Newton nonlinear solver class header file.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#ifndef NEWTON_DEFINED__
#define NEWTON_DEFINED__

// Inclusions
#include <math.h>
#include "matrix.hpp"
#include "resid.hpp"


// Newton solver class
class NewtonSolver {

 private:

  // private reusable data
  std::vector<double> f;      // stores nonlinear residual vector
  std::vector<double> s;      // stores Newton update vector
  Matrix J;                   // stores nonlinear residual Jacobian matrix

  // private pointers to problem-defining function objects
  ResidualFunction *fres;   // nonlinear residual function pointer
  ResidualJacobian *Jres;   // nonlinear residual Jacobian function pointer

  // private solver parameters
  double rtol;              // desired relative solution accuracy (non-negative)
  double atol;              // desired absolute solution accuracy (non-negative)
                            //   at most one of rtol,atol can equal zero
  int maxit;                // maximum desired Newton iterations
  bool show_iterates;       // flag for printing iteration information to screen

  // private statistics 
  int iters;                // iteration counter (reset in each solve)
  double error_norm;        // norm of most recent error estimate

 public:

  // Constructor
  NewtonSolver(ResidualFunction& fres_, ResidualJacobian& Jres_, 
	       const double rtol_, const double atol_, const int maxit_, 
               const std::vector<double>& y, const bool show_iterates_);

  // Newton solver routine
  int Solve(std::vector<double>& y);

  // Parameter update routines
  int SetMaxit(const int maxit_) {
    if (maxit_ < 1) { std::cerr << "SetMaxit error: illegal maxit\n";  return -1; }
    maxit = maxit_;
    return 0;
  };
  void SetShowIterates(const bool TF) { 
        show_iterates = TF; };
  int SetTolerances(const double rtol_, const double atol_) { 
    if ((rtol_ < 0.0) || (atol_ < 0.0) || (rtol_ == 0.0  && atol_ == 0.0)) {
      std::cerr << "SetTolerances error: illegal tolerances\n";
      return -1;
    }
    rtol = rtol_;
    atol = atol_;
    return 0;
  }
  void ResetIters() { iters = 0; };

  // Statistics accessor routines
  const int GetIters() { return iters; };
  const double GetErrorNorm() { return error_norm; };

};

#endif
