/* Linear multistep time stepper class header file. Unaltered from Dan
 * Reynolds' original file.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#ifndef LMM_DEFINED__
#define LMM_DEFINED__

// Inclusions
#include <vector>
#include <math.h>
#include "matrix.hpp"
#include "rhs.hpp"
#include "resid.hpp"
#include "newton.hpp"


// Linear multistep residual function class -- implements a 
// LMM-specific ResidualFunction to be supplied to the Newton solver.
class LMMResid: public ResidualFunction {
public:

  // data required to evaluate LMM nonlinear residual
  RHSFunction *frhs;       // pointer to ODE RHS function
  double t;                // current time
  double h;                // current step size
  Matrix yold;             // matrix of old solution vectors
  Matrix fold;             // matrix of old right-hand side vectors
  std::vector<double> a;   // vector of LMM "a" coefficients
  std::vector<double> b;   // vector of LMM "b" coefficients

  // constructor (sets RHS function and old solution vector pointers)
  LMMResid(RHSFunction& frhs_, std::vector<double> y, 
           std::vector<double> &a_, std::vector<double>& b_) {
    frhs = &frhs_;                      // set RHSFunction pointer
    a = a_;  b = b_;                    // copy LMM coefficients
    yold = Matrix(y.size(), a.size());  // allocate LMM arrays
    fold = Matrix(y.size(), a.size());
  };

  // initializer (fills initial set of 'old' vectors)
  int Initialize(double t, double h, Matrix& y);

  // residual evaluation routine
  int Evaluate(std::vector<double>& y, std::vector<double>& resid);

  // updater (shifts 'old' vectors, adding new one)
  int Update(double t, std::vector<double>& ynew);

};


// Linear multistep residual Jacobian function class -- implements 
// a LMM-specific ResidualJacobian to be supplied to the Newton solver.
class LMMResidJac: public ResidualJacobian {
public:

  // data required to evaluate LMM residual Jacobian
  RHSJacobian *Jrhs;   // ODE RHS Jacobian function pointer
  double t;            // current time
  double h;            // current step size
  double beta;         // b_{-1} coefficient

  // constructor (sets RHS Jacobian function pointer)
  LMMResidJac(RHSJacobian& Jrhs_, double beta_) { 
    Jrhs = &Jrhs_;   // set RHSJacobian pointer
    beta = beta_;    // copy b_{-1} coefficient
  };

  // Residual Jacobian evaluation routine
  int Evaluate(std::vector<double>& y, Matrix& J);
};



// LMM time stepper class
class LMMStepper {

 private:

  // private reusable local data
  RHSFunction *frhs;      // pointer to ODE RHS function
  LMMResid *resid;        // pointer to LMM residual function
  LMMResidJac *residJac;  // pointer to LMM residual Jacobian function
  std::vector<double> a;  // LMM coefficients
  std::vector<double> b;

 public:

  NewtonSolver *newt;     // Newton nonlinear solver pointer

  // constructor (constructs residual, Jacobian, and copies y for local data)
  LMMStepper(RHSFunction& frhs_, RHSJacobian& Jrhs_, std::vector<double>& y, 
             std::vector<double>& a_, std::vector<double>& b_);

  // destructor (frees local data)
  ~LMMStepper() {
    delete newt;
    delete resid;
    delete residJac;
  };

  // Evolve routine (evolves the solution via LMM)
  std::vector<double> Evolve(std::vector<double>& tspan, double h, Matrix& y);

};

#endif
