/* Adaptive forward Euler time stepper class header file.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#ifndef ADAPT_EULER_DEFINED__
#define ADAPT_EULER_DEFINED__

// Inclusions
#include <vector>
#include <math.h>
#include "matrix.hpp"
#include "rhs.hpp"


// Adaptive forward Euler time stepper class
class AdaptEuler {

 private:

  // private reusable local data
  RHSFunction *frhs;          // pointer to ODE RHS function
  std::vector<double> fn;     // local vector storage
  std::vector<double> y1;
  std::vector<double> y2;
  std::vector<double> yerr;

 public:

  double rtol;        // desired relative solution error
  double atol;        // desired absolute solution error
  double grow;        // maximum step size growth factor
  double safe;        // safety factor for step size estimate
  double fail;        // failed step reduction factor
  double ONEMSM;      // safety factors for
  double ONEPSM;      // floating-point comparisons
  double alpha;       // exponent relating step to error
  double error_norm;  // current estimate of the local error ratio
  double h;           // current time step size
  long int fails;     // number of failed steps
  long int steps;     // number of successful steps
  long int maxit;     // maximum number of steps

  // constructor (sets RHS function pointer & solver parameters, copies y for local data)
  AdaptEuler(RHSFunction& frhs_, double rtol_, double atol_, std::vector<double>& y);

  // Evolve routine (evolves the solution via adaptive forward Euler)
  std::vector<double> Evolve(std::vector<double>& tspan, std::vector<double>& y);

};

#endif
