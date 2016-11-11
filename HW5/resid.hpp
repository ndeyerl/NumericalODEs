/* Newton nonlinear solver class header file.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2014  */

#ifndef RESIDUAL_DEFINED__
#define RESIDUAL_DEFINED__

// Inclusions
#include <vector>
#include <math.h>
#include "matrix.hpp"


// Declare abstract base classes for residual and Jacobian to 
// define what the Newton solver expects from each.

//   Residual function abstract base class; derived classes 
//   must at least implement the Evaluate() routine
class ResidualFunction {
 public: 
  virtual int Evaluate(std::vector<double>& y, std::vector<double>& r) = 0;
};

//   Residual Jacobian function abstract base class; derived classes 
//   must at least implement the Evaluate() routine
class ResidualJacobian {
 public: 
  virtual int Evaluate(std::vector<double>& y, Matrix& J) = 0;
};

#endif
