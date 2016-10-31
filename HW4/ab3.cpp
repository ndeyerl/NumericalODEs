/* Explicit 4th-order Adams-Bashforth time stepper class implementation file.
 * Based on Dan Reynolds ab1.cpp file.

   Nicole Deyerl
   Math 6321 @ SMU
   Fall 2016  */

#include <vector>
#include "matrix.hpp"
#include "ab3.hpp"


// The explicit 4th-order Adams-Bashforth time step evolution routine
//
// Inputs:  tspan holds the current time interval, [t0, tf]
//          h holds the desired time step size
//          y holds the initial condition, y(t0) or y(n)
//          y1 holds the previous initial condition, y(t0-h) or y(n-1)
//          y2 holds y(t0-2*h) or y(n-2)
//          y3 holds y(t0-3*h) or y(n-3)
// Outputs: y holds the computed solution, y(tf)
//          y1 holds the next-to-last computed solution, y(tf-h)
//          y2 holds the computed solution y(tf-2*h)
//          y3 holds the computed solution y(tf-3*h)
//
// The return value is a row vector containing all internal 
// times at which the solution was computed,
//               [t0, t1, ..., tN]
std::vector<double> AB3Stepper::Evolve(std::vector<double>& tspan, 
			     double h, 
                 std::vector<double>& y, 
			     std::vector<double>& y1,
			     std::vector<double>& y2,
			     std::vector<double>& y3) {

  // initialize output
  std::vector<double> times = {tspan[0]};

  // check for legal inputs 
  if (h <= 0.0) {
    std::cerr << "Evolve: Illegal h\n";
    return times;
  }
  if (tspan[1] <= tspan[0]) {
    std::cerr << "Evolve: Illegal tspan\n";
    return times;	  
  }
  
  // figure out how many time steps
  long int N = (tspan[1]-tspan[0])/h;
  if (tspan[1] > tspan[0]+N*h)  N++;

  //evaluate frhs at y1, y2, y3 and store in f, f1, f2
  // store these shifted over so that the update in the 
  // iteration is preserved and works for first step and 
  // all steps afterwards
  if (frhs->Evaluate(tspan[0]-h, y1, f) != 0) {
    std::cerr << "Evolve: Error in ODE RHS function\n";
    return times;
  }
  if (frhs->Evaluate(tspan[0]-2.0*h, y2, f1) != 0) {
    std::cerr << "Evolve: Error in ODE RHS function\n";
    return times;
  }
  if (frhs->Evaluate(tspan[0]-3.0*h, y3, f2) != 0) {
    std::cerr << "Evolve: Error in ODE RHS function\n";
    return times;
  }    
    
  // iterate over time steps
  for (int i=0; i<N; i++) {

    // last step only: update h to stop directly at final time
    // NOTE: if  this actually differs from the input h, then 
    //       the LMM will reduce to 1st order
    if (i == N-1) 
      h = tspan[1]-times[i];

    // update old f's and y's
    y3 = y2;
    y2 = y1;
    y1 = y;
    f3 = f2;
    f2 = f1;
    f1 = f;

    // evaluate f at the current y
    if (frhs->Evaluate(times[i], y, f) != 0) {
      std::cerr << "Evolve: Error in ODE RHS function\n";
      return times;
    }

    // update the current solution using AB3 method from prob 1
    y += (h/24.0)*(55.0*f-59.0*f1+37.0*f2-9.0*f3);

    // update current time, store in output array
    times.push_back(times[i] + h);
  }

  return times;
}
