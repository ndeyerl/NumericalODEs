/* Adaptive explicit RKF-45 time stepper class header file.
 * Based off of Dan Reynolds' adapt_euler and RK4 files.

   Nicole Deyerl
   Math 6321 @ SMU
   Fall 2016  */

#ifndef ADAPT_RKF_DEFINED__
#define ADAPT_RKF_DEFINED__

// Inclusions
#include <vector>
#include <math.h>
#include "matrix.hpp"
#include "rhs.hpp"


// Adaptive explicit RKF-45 time stepper class
class AdaptRKF {

 private:

  // private reusable local data
  RHSFunction *frhs;          // pointer to ODE RHS function
  std::vector<double> fn;     // local vector storage
  std::vector<double> y4;
  std::vector<double> y5;
  std::vector<double> yerr;
  std::vector<double> z, f0, f1, f2, f3, f4, f5;    // reused RK vectors
  Matrix A;                                 // Butcher table
  std::vector<double> b5, b4, c; //RK5 and RK4 coeffs

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
  AdaptRKF(RHSFunction& frhs_, double rtol_, double atol_, std::vector<double>& y) {
    frhs = &frhs_;                        // store RHSFunction pointer
    z  = y;   // allocate reusable data 
    f0 = y;   //   based on size of y
    f1 = y;
    f2 = y;
    f3 = y;
    f4 = y;
    f5 = y;
    A = Matrix(6,6);                      // Butcher table data
    A(1,0) = 0.25;
    A(2,0) = 0.09375; A(2,1) = 0.28125;
    A(3,0) = 1932.0/2197.0;  A(3,1) = -7200.0/2197.0; A(3,2) = 7296.0/2197.0;
    A(4,0) = 439.0/216.0;  A(4,1) = -8.0;  A(4,2) = 3680.0/513.0;  A(4,3) = -845.0/4104.0;
    A(5,0) = -8.0/27.0;  A(5,1) = 2.0;  A(5,2) = -3544.0/2565.0;  A(5,3) = 1859.0/4104.0; A(5,4) = -11.0/40.0;
    b5 = {16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0};
    b4 = {25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0}; 
    c = {0.0, 0.25, 0.375, 12.0/13.0, 1.0, 0.5};
    
    rtol = rtol_;     // set tolerances
    atol = atol_;
    fn = y;           // clone y to create local vectors
    y4 = y;
    y5 = y;
    yerr = y;

    maxit = 1e6;      // set default solver parameters
    grow = 50.0;
    safe = 0.95;
    fail = 0.5;
    ONEMSM = 1.0 - 1.e-8;
    ONEPSM = 1.0 + 1.e-8;
    alpha = -0.5;
    fails = 0;
    steps = 0;
    error_norm = 0.0;
    h = 0.0;
  };    

  // Evolve routine (evolves the solution via adaptive RKF45)
  std::vector<double> Evolve(std::vector<double>& tspan, std::vector<double>& y);

  // Single RK4, RK5 step calculation
  int Step(double t, double h, std::vector<double>& y, std::vector<double>& y4, 
                              std::vector<double>& y5);
  
};

#endif
