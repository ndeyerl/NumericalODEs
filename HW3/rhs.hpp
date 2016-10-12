/* ODE RHS and Jacobian function abstract base class definitions.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#ifndef ODE_RHS_DEFINED__
#define ODE_RHS_DEFINED__

// Inclusions
#include "matrix.hpp"


// Declare abstract base classes for ODE RHS and its Jacobian, to 
// define what the backward Euler solver expects from each.

//   ODE RHS function abstract base class; derived classes 
//   must at least implement the Evaluate() routine
class RHSFunction {
 public: 
  virtual int Evaluate(double t, std::vector<double>& y, std::vector<double>& f) = 0;
};

//   ODE RHS Jacobian function abstract base class; derived 
//   classes must at least implement the Evaluate() routine
class RHSJacobian {
 public: 
  virtual int Evaluate(double t, std::vector<double>& y, Matrix& J) = 0;
};

#endif
