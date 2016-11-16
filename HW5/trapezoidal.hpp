/* Trapezoidal time stepper class header file.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#ifndef TRAPEZOIDAL_DEFINED__
#define TRAPEZOIDAL_DEFINED__

// Inclusions
#include <vector>
#include <math.h>
#include "matrix.hpp"
#include "rhs.hpp"
#include "resid.hpp"
#include "newton.hpp"


// Backward Euler residual function class -- implements a 
// backward-Euler-specific ResidualFunction to be supplied 
// to the Newton solver.
class TrapResid: public ResidualFunction {
public:

  // data required to evaluate backward Euler nonlinear residual
  RHSFunction *frhs;            // pointer to ODE RHS function
  double t;                     // current time
  double h;                     // current step size
  std::vector<double> *yold;    // pointer to solution at old time step
  std::vector<double>  fold;    // extra vector for residual evaluation

  // constructor (sets RHSFunction and old solution pointers)
  TrapResid(RHSFunction& frhs_, std::vector<double>& yold_) {
    frhs = &frhs_;  yold = &yold_;  fold = yold_;
  };

  // residual evaluation routine
  int Evaluate(std::vector<double>& y, std::vector<double>& resid) {

    // evaluate RHS function at new time (store in resid)
    int ierr = frhs->Evaluate(t+h, y, resid);
    if (ierr != 0) {
      std::cerr << "Error in ODE RHS function = " << ierr << "\n";
      return ierr;
    }

    // evaluate RHS function at old time
    ierr = frhs->Evaluate(t, (*yold), fold);
    if (ierr != 0) {
      std::cerr << "Error in ODE RHS function = " << ierr << "\n";
      return ierr;
    }

    // combine pieces to fill residual, y-yold-h/2*[f(t+h,y)+f(t,yold)]
    resid = y - (*yold) - 0.5*h*(resid+fold);

    // return success
    return 0;
  }
};



// Backward Euler residual Jacobian function class -- implements 
// a backward-Euler-specific ResidualJacobian to be supplied 
// to the Newton solver.
class TrapResidJac: public ResidualJacobian {
public:

  // data required to evaluate backward Euler residual Jacobian
  RHSJacobian *Jrhs;   // ODE RHS Jacobian function pointer
  double t;            // current time
  double h;            // current step size

  // constructor (sets RHS Jacobian function pointer)
  TrapResidJac(RHSJacobian &Jrhs_) { Jrhs = &Jrhs_; };

  // Residual Jacobian evaluation routine
  int Evaluate(std::vector<double>& y, Matrix& J) {

    // evaluate RHS function Jacobian, Jrhs (store in J)
    int ierr = Jrhs->Evaluate(t+h, y, J);
    if (ierr != 0) {
      std::cerr << "Error in ODE RHS Jacobian function = " << ierr << "\n";
      return ierr;
    }
    // combine pieces to fill residual Jacobian,  J = I - h/2*Jrhs
    J *= (-0.5*h);
    for (int i=0; i<J.Rows(); i++)
      J(i,i) += 1.0;

    // return success
    return 0;
  }
};



// Trapezoidal time stepper class
class TrapezoidalStepper {

 private:

  // private reusable local data
  std::vector<double> yold;   // old solution vector
  TrapResid *resid;           // trapezoidal Euler residual function pointer
  TrapResidJac *residJac;     // trapezoidal residual Jacobian function pointer

 public:

  NewtonSolver *newt;   // Newton nonlinear solver pointer

  // constructor (constructs residual, Jacobian, and copies y for local data)
  TrapezoidalStepper(RHSFunction& frhs_, RHSJacobian& Jrhs_, 
                     std::vector<double>& y_);

  // destructor (frees pointers to local objects)
  ~TrapezoidalStepper() {
    delete resid;
    delete residJac;
    delete newt;
  };

  // Evolve routine (evolves the solution via trapezoidal)
  std::vector<double> Evolve(std::vector<double>& tspan, double h, 
                             std::vector<double>& y);

};

#endif
