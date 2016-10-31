/* Adaptive Runge Kutta Fehlberg (45) solver class implementation file.
   Based off of Dan Reynolds' adapt_euler and RK45 scripts.

   Class to perform time evolution of the IVP
        y' = f(t,y),  t in [t0, Tf],  y(t0) = y0
   using the RKF-45 time stepping method. 

   Nicole Deyerl
   Math 6321 @ SMU
   Fall 2016  */

#include "matrix.hpp"
#include "adapt_rkf.hpp"


// The adaptive RKF45 time step evolution routine
//
// Inputs:  tspan holds the current time interval, [t0, tf]
//          y holds the initial condition, y(t0)
// Outputs: y holds the computed solution, y(tf)
//
// The return value is a row vector containing all internal 
// times at which the solution was computed,
//               [t0, t1, ..., tN]
std::vector<double> AdaptRKF::Evolve(std::vector<double>& tspan, 
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
    y4 = y;
    y5 = y;

    // perform a single step of RK45 to update y
    if (Step(tvals[steps], h, y, y4, y5) != 0) {
      std::cerr << "Evolve: Error in Step() function\n";
      return tvals;
    }
    
    // compute error estimate
    yerr = y5 - y4;

    // compute error estimate success factor
    error_norm = std::max(InfNorm(yerr) / ( rtol * InfNorm(y5) + atol ), 1.e-8);

    // if solution has too much error: reduce step size, increment failure counter, and retry
    if (error_norm > ONEPSM) {
      h *= fail;
      fails++;
      continue;
    }

    // successful step
    tvals.push_back(tvals[steps++] + h);  // append updated time, increment step counter
    y = y4;  //in general we keep the "worse" step, not the better one; the
             // better one is used solely as the error estimator

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

// Single step of explicit 4th and 5th-order Runge-Kutta
//
// Inputs:  t holds the current time
//          h holds the current time step size
//          z, f0-f5 hold temporary vectors needed for the problem
//          y holds the current solution (from RK4)
// Outputs: y4 holds the updated solution, y(t+h) for RK4
//          y5 holds the updated solution, y(t+h) for RK5
//
// The return value is an integer indicating success/failure,
// with 0 indicating success, and nonzero failure.
int AdaptRKF::Step(double t, double h, std::vector<double>& y,
                      std::vector<double>& y4, std::vector<double>& y5) {

  // stage 1: set stage and compute RHS
  z = y;
  if (frhs->Evaluate(t, z, f0) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  }

  // stage 2: set stage and compute RHS
  z = y + (h*A(1,0))*f0;
  if (frhs->Evaluate(t+c[1]*h, z, f1) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  }

  // stage 3: set stage and compute RHS
  z = y + h*(A(2,0)*f0 + A(2,1)*f1);
  if (frhs->Evaluate(t+c[2]*h, z, f2) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  }

  // stage 4: set stage and compute RHS
  z = y + h*(A(3,0)*f0 + A(3,1)*f1 + A(3,2)*f2);
  if (frhs->Evaluate(t+c[3]*h, z, f3) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  }
  
  // stage 5: set stage and compute RHS
  z = y + h*(A(4,0)*f0 + A(4,1)*f1 + A(4,2)*f2 + A(4,3)*f3);
  if (frhs->Evaluate(t+c[4]*h, z, f4) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  }
  
  // stage 6: set stage and compute RHS
  z = y + h*(A(5,0)*f0 + A(5,1)*f1 + A(5,2)*f2 + A(5,3)*f3 + A(5,4)*f4);
  if (frhs->Evaluate(t+c[5]*h, z, f5) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  }

  // compute next step solutions
  y4 += h*(b4[0]*f0 + b4[1]*f1 + b4[2]*f2 + b4[3]*f3 + b4[4]*f4);
  y5 += h*(b5[0]*f0 + b5[1]*f1 + b5[2]*f2 + b5[3]*f3 + b5[4]*f4 + b5[5]*f5);

  // return success
  return 0;
};
