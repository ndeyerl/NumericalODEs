/* Main routine to test the fourth-order Adams 
   Bashforth method on the scalar-valued ODE problem 
     y' = -exp(-t)*y, t in [0,5],
     y(0) = 1.
   
   Based off of Dan Reynolds ERK driver.

   Nicole Deyerl
   Math 6321 @ SMU
   Fall 2016  */

#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"
#include "erk4.hpp"
#include "ab3.hpp"

using namespace std;


// ODE RHS function
class MyRHS: public RHSFunction {
public:
  int Evaluate(double t, vector<double>& y, vector<double>& f) {
    f[0] = -exp(-t)*y[0];
    return 0;
  }
};


// Convenience function for analytical solution
vector<double> ytrue(const double t) { 
  vector<double> yt = {exp(exp(-t)-1.0)};
  return yt;
};


// main routine
int main() {

  // time steps to try
  vector<double> h = {0.1, 0.05, 0.01, 0.005, 0.001, 0.0005};

  // set problem information
  vector<double> y0 = {1.0};
  double t0 = 0.0;
  double Tf = 5.0;
  double dtout = 1.0;

  // create ODE RHS function objects
  MyRHS rhs;

  // create ERK4 and AB3 stepper objects
  ERK4Stepper ERK4(rhs, y0);
  AB3Stepper AB3(rhs, y0);

  // storage for error results
  vector<double> errs(h.size());

  ///////// Adams Bashforth 3 /////////
  cout << "\nAdams Bashforth 3 Method:\n";

  // loop over time step sizes
  for (int ih=0; ih<h.size(); ih++) {

    cout << "Running problem 2 with h = " << h[ih] << endl;
    
    // set the initial conditions (using RK4) and initial time
    vector<double> y(y0);
    double tcur = t0;
    vector<double> tspan = {tcur,tcur+h[ih]};
    vector<double> tvals = ERK4.Evolve(tspan, h[ih], y);
    tcur = tvals.back();
    vector<double> y3 = y;
    tspan = {tcur,tcur+h[ih]};
    tvals = ERK4.Evolve(tspan, h[ih], y);
    tcur = tvals.back();
    vector<double> y2 = y;
    tspan = {tcur,tcur+h[ih]};
    tvals = ERK4.Evolve(tspan, h[ih], y);
    tcur = tvals.back();
    vector<double> y1 = y;
    tspan = {tcur, tcur+h[ih]};
    tvals = ERK4.Evolve(tspan, h[ih], y);
    tcur = tvals.back();

    // reset maxerr
    double maxerr = 0.0;
    
    // variables to handle the timing being off from RK4 IC's
    double tf;
    int count = 0;
    // loop over output step sizes: call solver and output error
    while (tcur < 0.99999*Tf) {
      
      // set the time interval for this solve, fix it for the first
      // dtout
      if ( count == 0 ){
          tf = tcur-4.0*h[ih] + dtout;
      } else {
		  tf = tcur + dtout;
	  }
	  count = count + 1;
      vector<double> tspan = {tcur, std::min(tf, Tf)};

      // call the solver, update current time
      vector<double> tvals = AB3.Evolve(tspan, h[ih], y, y1, y2, y3);
      tcur = tvals.back();     // last entry in tvals
      
      // compute the abs error at tcur, output to screen and accumulate maximum
      vector<double> yerr = y - ytrue(tcur);
      double err = InfNorm(yerr);
      maxerr = std::max(maxerr, err); // holds max abs error
      cout << " y(" << tcur << ") = " << y[0]
      << " abserr = " << err << endl;
      
    }
    cout << "   h = " << h[ih] << "\t  maxerror = " << maxerr;
    errs[ih] = maxerr;
    
    //convergence print out
    if (ih > 0)
      cout << "\t  conv rate = " << (log(errs[ih])-log(errs[ih-1]))/(log(h[ih])-log(h[ih-1]));
    cout << endl;
  }


  return 0;
}
