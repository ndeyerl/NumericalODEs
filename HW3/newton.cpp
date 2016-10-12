/* Newton nonlinear solver class implementation file.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#include "matrix.hpp"
#include "newton.hpp"
#include "resid.hpp"


// Newton solver construction routine
//
// Inputs:  fres_  -- the ResidualFunction to use
//          Jres_  -- the JacobianFunction to use
//          rtol_  -- the desired solution residual tolerance (non-negative)
//          atol_  -- the desired solution absolute tolerance (non-negative)
//          maxit_ -- the maximum allowed number of iterations
//          y      -- an example solution vector (only used for size/shape)
NewtonSolver::NewtonSolver(ResidualFunction& fres_, ResidualJacobian& Jres_, 
                           const double rtol_, const double atol_, const int maxit_, 
                           const std::vector<double>& y, const bool show_iterates_) {

  // set pointers to problem-defining function objects
  fres = &fres_;
  Jres = &Jres_;

  // set solver parameters
  SetShowIterates(show_iterates_);
  if (SetMaxit(maxit_) != 0) {
    std::cerr << "cannot create NewtonSolver object"; 
    return;
  }
  if (SetTolerances(rtol_, atol_) != 0) { 
    std::cerr << "cannot create NewtonSolver object"; 
    return;
  }

  // create reusable solver objects (clone off of y)
  f = std::vector<double>(y.size());
  s = std::vector<double>(y.size());
  J = Matrix(y.size(), y.size());
  
  // initialize statistics
  iters = 0;
  error_norm = 0.0;

};



// The actual Newton solver routine
//
// Input:   y  -- the initial guess
// Outputs: y  -- the computed solution
//  
// The return value is one of:
//          0 => successful solve
//         -1 => bad function call or input
//          1 => non-convergent iteration
int NewtonSolver::Solve(std::vector<double>& y) {

  // set initial residual value
  if (fres->Evaluate(y, f) != 0) {
    std::cerr << "NewtonSolver::Solve error: residual function failure\n";
    return -1;
  }


  // perform iterations
  for (iters=1; iters<=maxit; iters++) {

    // evaluate Jacobian 
    if (Jres->Evaluate(y, J) != 0) {
      std::cerr << "NewtonSolver::Solve error: Jacobian function failure\n";
      return -1;
    }

    // compute Newton update, norm
    if (LinearSolve(J, s, f) != 0) {
      std::cerr << "NewtonSolver::Solve error: linear solver failure\n";
      return -1;
    }
    error_norm = InfNorm(s);

    // perform update
    y -= s;

    // update residual
    if (fres->Evaluate(y, f) != 0) {
      std::cerr << "NewtonSolver::Solve error: residual function failure\n";
      return -1;
    }

    // output convergence information
    if (show_iterates)
      std::cout << "   iter " << iters
                << ", ||s||_inf = " << error_norm
                << ", ||f(x)||_inf = " << InfNorm(f) << std::endl;

    // check for convergence, return if successful
    if (error_norm < rtol*InfNorm(y) + atol) 
    return 0;

  }

  // if we've made it here, Newton did not converge, output warning and return
  std::cerr << "\nNewtonSolver::Solve WARNING: nonconvergence after " << maxit 
	    << " iterations (||s|| = " << error_norm << ")\n";
  return 1;
}

