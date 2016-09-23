/* Forward Euler time stepper class header file.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#ifndef FORWARD_EULER_DEFINED__
#define FORWARD_EULER_DEFINED__

// Inclusions
#include <vector>
#include <math.h>
#include "matrix.hpp"
#include "rhs.hpp"


// Forward Euler time stepper class
class ForwardEulerStepper {

 private:

  // private reusable local data
  std::vector<double> f;   // storage for ODE RHS vector
  RHSFunction *frhs;       // pointer to ODE RHS function

 public:

  // constructor (sets RHS function pointer, copies y for local data)
  ForwardEulerStepper(RHSFunction& frhs_, std::vector<double>& y) {
    frhs = &frhs_;
    f = y;
  };

  // Evolve routine (evolves the solution via forward Euler)
  std::vector<double> Evolve(std::vector<double> tspan, double h, 
                             std::vector<double>& y);

};

#endif
