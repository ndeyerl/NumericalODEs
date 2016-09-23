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
class AdaptEuler {

 private:

  // private reusable local data
  std::vector<double> f;   // storage for ODE RHS vector
  RHSFunction *frhs;       // pointer to ODE RHS function
  double r; //storage for rtol
  double a; //storage for atol

 public:

  // constructor (sets RHS function pointer, copies y for local data)
  AdaptEuler(RHSFunction& frhs_, double rtol, double atol, std::vector<double>& y) {
    frhs = &frhs_;
    r = rtol;
    a = atol;
    f = y;
  };

  // Evolve routine (evolves the solution via forward Euler)
  std::vector<double> Evolve(std::vector<double> tspan, std::vector<double>& y);

};

#endif
