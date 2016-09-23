/* Forward Euler time stepper class implementation file.

   Class to perform time evolution of the IVP
        y' = f(t,y),  t in [t0, Tf],  y(t0) = y0
   using the forward Euler (explicit Euler) time stepping method. 

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016 */

#include <vector>
#include "matrix.hpp"
#include "fwd_euler.hpp"


// The forward Euler time step evolution routine
//
// Inputs:  tspan holds the current time interval, [t0, tf]
//          h holds the desired time step size
//          y holds the initial condition, y(t0)
// Outputs: y holds the computed solution, y(tf)
//
// The return value is a row vector containing all internal 
// times at which the solution was computed,
//               [t0, t1, ..., tN]
std::vector<double> ForwardEulerStepper::Evolve(std::vector<double> tspan, double h, std::vector<double>& y) {

  // initialize output
  std::vector<double> times = {tspan[0]};

  // check for legal inputs 
  if (h <= 0.0) {
    std::cerr << "ForwardEulerStepper: Illegal h\n";
    return times;
  }
  if (tspan[1] <= tspan[0]) {
    std::cerr << "ForwardEulerStepper: Illegal tspan\n";
    return times;	  
  }
  
  // figure out how many time steps
  long int N = (tspan[1]-tspan[0]) / h;
  if (tspan[1] > tspan[0]+N*h)  N++;
    
  // iterate over time steps
  for (long int i=0; i<N; i++) {

    // last step only: update h to stop directly at final time
    if (i == N-1) 
      h = tspan[1]-times[i];

    // compute ODE RHS
    if (frhs->Evaluate(times[i], y, f) != 0) {
      std::cerr << "ForwardEulerStepper: Error in ODE RHS function\n";
      return times;
    }

    // update solution with forward Euler step
    y += (h*f);

    // update current time, store in output array
    times.push_back(times[i] + h);
  }

  return times;
}
