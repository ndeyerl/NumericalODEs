/* Generalized Trapezoidal time stepper class implementation file.

   Class to perform time evolution of the IVP
        y' = f(t,y),  t in [t0, Tf],  y(t0) = y0
   using the generalized trapezoidal time stepping method. 

   Nicole Deyerl
   Math 6321 @ SMU
   Fall 2016  */ 

#include "trapezoidal.hpp" 


// Trapezoidal stepper construction routine (allocates local data)
//
// Inputs:  frhs_  holds the RHSFunction to use
//          Jrhs_  holds the RHSJacobian to use
//          theta  holds the value for theta for desired trapezoid method
//          y      holds an example solution vector (only used for cloning)
TrapezoidalStepper::TrapezoidalStepper(RHSFunction& frhs_, RHSJacobian& Jrhs_, double theta, 
                                       std::vector<double>& y) {

  // clone y to create local reusable data
  yold = y; 

  // construct objects for nonlinear residual and its Jacobian
  resid = new TrapResid(frhs_, yold);
  residJac = new TrapResidJac(Jrhs_);
  
  // copy current theta value into resid and resid jac objects
  resid->theta = theta;
  residJac->theta = theta;

  // construct Newton solver object (only copies y)
  // (initialize with default solver parameters; user may override with, e.g.
  //  TrapezoidalStepper.newt->SetMaxit())
  newt = new NewtonSolver(*resid, *residJac, 1.0e-7, 1.0e-11, 100, y, false);
};


// The actual trapezoidal time step evolution routine
//
// Inputs:  tspan holds the current time interval, [t0, tf]
//          h holds the desired time step size
//          y holds the initial condition, y(t0)
// Outputs: y holds the computed solution, y(tf)
//
// The return value is a row vector containing all internal 
// times at which the solution was computed,
//               [t0, t1, ..., tN]
std::vector<double> TrapezoidalStepper::Evolve(std::vector<double>& tspan, double h, 
                                               std::vector<double>& y) {

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

  // iterate over time steps
  for (int i=0; i<N; i++) {

    // last step only: update h to stop directly at final time
    if (i == N-1) 
      h = tspan[1]-times[i];

    // update resid and residJac objects with information on current state
    resid->t    = times[i];    // copy current time into objects
    residJac->t = times[i];
    resid->h    = h;           // copy current stepsize into objects
    residJac->h = h;
    yold = y;                  // copy y into stored yold object

    // call Newton method to solve for the updated solution
    int ierr = newt->Solve(y);
    if (ierr != 0) {
      std::cerr << "TrapezoidalStepper: Error in Newton solver function = " 
                << ierr << "\n";
      return times;
    }

    // update current time, store in output array
    times.push_back(times[i] + h);
  }

  return times;
}
