/* This is a routine to test the forward Euler method for a scalar valued ODE.
     y' = f(t,y), t in [0,5],
     y(0) = y0.
   It is based off of Dan Reynold's driver_fwd_euler and uses his suite of 
   forward Euler solvers.  


   Nicole Deyerl
   Math 6321 
   Fall 2016  */

#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"
#include "fwd_euler.hpp"

using namespace std;


// Problem 1:  y' = -y
//    ODE RHS function class -- instantiates a RHSFunction
class RHS1: public RHSFunction {
public:
  int Evaluate(double t, vector<double>& y, vector<double>& f) {    // evaluates the RHS function, f(t,y)
    f = -exp(-t)*y;
    return 0;
  }
};
// Convenience function for analytical solution
//   y(t) = e^(e^(-t)-1)
vector<double> ytrue1(const double t) { 
  vector<double> yt(1);
  yt[0] = exp(exp(-t)-1);
  return yt;
};



// Problem 2:  y' = (y+t^2-2)/(t+1)
//    ODE RHS function class -- instantiates a RHSFunction
class RHS2: public RHSFunction {
public:
  int Evaluate(double t, vector<double>& y, vector<double>& f) {    // evaluates the RHS function, f(t,y)
    f[0] = (y[0]+t*t-2.0)/(t+1.0);
    return 0;
  }
};
//    Convenience function for analytical solution
vector<double> ytrue2(const double t) { 
  vector<double> yt(1);
  yt[0] = t*t + 2.0*t + 2.0 - 2.0*(t+1.0)*log(t+1.0);
  return yt;
};



// main routine
int main() {

  // time steps to try
  vector<double> h = {0.2,0.1,0.05,0.025,0.0125}; //specified by 
                                                  // problem 1a
  // set problem information
  vector<double> y0_1 = {1.0};
  vector<double> y0_2 = {2.0};
  double t0 = 0.0;
  double Tf = 5.0;
  double tcur = t0;
  double dtout = 1.0;

  // create ODE RHS function objects
  RHS1 f1;
  RHS2 f2;

  // create forward Euler stepper object
  ForwardEulerStepper FE1(f1, y0_1);
  ForwardEulerStepper FE2(f2, y0_2);
  vector<double> maxerr1 (h.size(),0.0);
  vector<double> maxerr2 (h.size(),0.0);
  
  // loop over time step sizes
  for (int ih=0; ih<h.size(); ih++) {

    // problem 1:
    vector<double> y = y0_1;
    tcur = t0;
    cout << "\nRunning problem 1 with stepsize h = " << h[ih] << ":\n";

    //   loop over output step sizes: call solver and output error
    while (tcur < 0.99999*Tf) {
      
      // set the time interval for this solve
      vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};

      // call the solver, update current time
      vector<double> tvals = FE1.Evolve(tspan, h[ih], y);
      tcur = tvals.back();   // last entry in tvals

      // compute the error at tcur, output to screen and accumulate maximum
      vector<double> yerr = y - ytrue1(tcur);
      double err = InfNorm(yerr);
      maxerr1[ih] = std::max(maxerr1[ih], err);
      cout << "  y(" << tcur << ") = " << y[0]
	   << "  \t||error|| = " << err
	   << endl;
    }
    cout << "Max error = " << maxerr1[ih] << endl;


    // problem 2:
    y = y0_2;
    tcur = t0;
    cout << "\nRunning problem 2 with stepsize h = " << h[ih] << ":\n";

    //   loop over output step sizes: call solver and output error
    while (tcur < 0.99999*Tf) {
      
      // set the time interval for this solve
      vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};

      // call the solver, update current time
      vector<double> tvals = FE2.Evolve(tspan, h[ih], y);
      tcur = tvals.back();   // last entry in tvals

      // compute the error at tcur, output to screen and accumulate maximum
      vector<double> yerr = y - ytrue2(tcur);
      double err = InfNorm(yerr);
      maxerr2[ih] = std::max(maxerr2[ih], err);
      cout << "  y(" << tcur << ") = " << y[0]
	   << "  \t||error|| = " << err
	   << endl;
    }
    cout << "Max error = " << maxerr2[ih] << endl;

  }

  return 0;
}
