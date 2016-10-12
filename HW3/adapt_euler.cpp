/* Adaptive forward Euler solver class implementation file.

   Class to perform time evolution of the IVP
        y' = f(t,y),  t in [t0, Tf],  y(t0) = y0
   using the forward Euler (explicit Euler) time stepping method. 

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#include "matrix.hpp"
#include "adapt_euler.hpp"



// Adaptive forward Euler class constructor routine
//
// Inputs:  frhs_ holds the ODE RHSFunction object, f(t,y)
//          rtol holds the desired relative solution accuracy
//          atol holds the desired absolute solution accuracy
// 
// Sets default values for adaptivity parameters, all of which may
// be modified by the user after the solver object has been created
AdaptEuler::AdaptEuler(RHSFunction& frhs_, double rtol_, 
                       double atol_, std::vector<double>& y) {
  frhs = &frhs_;    // set RHSFunction pointer
  rtol = rtol_;     // set tolerances
  atol = atol_;
  fn = y;           // clone y to create local vectors
  y1 = y;
  y2 = y;
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


// The adaptive forward Euler time step evolution routine
//
// Inputs:  tspan holds the current time interval, [t0, tf]
//          y holds the initial condition, y(t0)
// Outputs: y holds the computed solution, y(tf)
//
// The return value is a row vector containing all internal 
// times at which the solution was computed,
//               [t0, t1, ..., tN]
std::vector<double> AdaptEuler::Evolve(std::vector<double>& tspan, 
                                       std::vector<double>& y) {

  // initialize output
  std::vector<double> tvals = {tspan[0]};

  // check for positive tolerances
  if ((rtol <= 0.0) || (atol <= 0.0)) {
    std::cerr << "Evolve error: illegal tolerances, atol = " 
	      << atol << ",  rtol = " << rtol << std::endl;
    return tvals;
  }

  // reset solver statistics
  fails = 0;
  steps = 0;

  // get |y'(t0)|
  if (frhs->Evaluate(tspan[0], y, fn) != 0) {
    std::cerr << "Evolve error in RHS function\n";
    return tvals;
  }

  // estimate initial h value via linearization, safety factor
  error_norm = std::max(Norm(fn) / ( rtol * Norm(y) + atol ), 1.e-8);
  h = safe/error_norm;
  if (tspan[0]+h > tspan[1])   h = tspan[1]-tspan[0];

  // iterate over time steps (all but the last one)
  for (int tstep=1; tstep<=maxit; tstep++) {

    // reset both solution approximations to current solution
    y1 = y;
    y2 = y;

    // get RHS at this time, perform full/half step updates
    if (frhs->Evaluate(tvals[steps], y, fn) != 0) {
      std::cerr << "Evolve error in RHS function\n";
      return tvals;
    }
    y1 += h*fn;
    y2 += (0.5*h)*fn;

    // get RHS at half-step, perform half step update
    if (frhs->Evaluate(tvals[steps]+0.5*h, y2, fn) != 0) {
      std::cerr << "Evolve error in RHS function\n";
      return tvals;
    }
    y2 += (0.5*h)*fn;

    // compute error estimate
    yerr = y2 - y1;

    // compute error estimate success factor
    error_norm = std::max(InfNorm(yerr) / ( rtol * InfNorm(y2) + atol ), 1.e-8);

    // if solution has too much error: reduce step size, increment failure counter, and retry
    if (error_norm > ONEPSM) {
      h *= fail;
      fails++;
      continue;
    }

    // successful step
    tvals.push_back(tvals[steps++] + h);  // append updated time, increment step counter
    y = 2.0*y2 - y1;

    // exit loop if we've reached the final time
    if ( tvals[steps] >= tspan[1]*ONEMSM )  break;

    // pick next time step size based on this error estimate
    double eta = safe*std::pow(error_norm, alpha);    // step size estimate
    eta = std::min(eta, grow);                        // maximum growth
    h *= eta;                                         // update h
    h = std::min(h, tspan[1]-tvals[steps]);           // don't pass Tf

  }

  // set output array as the subset of tvals that we actually used
  return tvals;
}
