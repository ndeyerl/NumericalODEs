/* Homework 4, problem 4: test adaptive RKF-45 and compare with adaptive 
 * Forward Euler for the scalar-valued ODE problem 
     y' = -y + 2cos(t), t in [0,10],
     y(0) = 1.
   Based off of Dan Reynolds' driver for hw2.

   Nicole Deyerl
   Math 6321 @ SMU
   Fall 2016  */

#include <iostream>
#include <iomanip>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"
#include "adapt_rkf.hpp"
#include "adapt_euler.hpp"

using namespace std;


// ODE RHS function class -- instantiates a RHSFunction
class MyRHS: public RHSFunction {
public:
  int Evaluate(double t, vector<double>& y, vector<double>& f) {
    f[0] = -y[0] + 2.0*cos(t);
    return 0;
  }
};
//    Convenience function for analytical solution
vector<double> ytrue(const double t) { 
  vector<double> yt(1);
  yt[0] = sin(t) + cos(t);
  return yt;
};


// main routine
int main() {

  // tolerances to try
  vector<double> rtols = {1.e-4, 1.e-6, 1.e-8};
  vector<double> atols = {1.e-6, 1.e-8, 1.e-10};

  // initial condition and time span
  vector<double> y0 = {1.0};
  double t0 = 0.0;
  double Tf = 10.0;
  double tcur = t0;
  double dtout = 1.0;

  // create ODE RHS function objects
  MyRHS rhs;

  // loop over relative tolerances
  for (int ir=0; ir<rtols.size(); ir++) {
 
    // create adaptive RKF-45 stepper object
    AdaptRKF ARKF(rhs, 0.0, atols[ir], y0);
    // set up the problem for this tolerance
    ARKF.rtol = rtols[ir];
    vector<double> y = y0;
    tcur = t0;
    double maxrelerr = 0.0;
    long int totsteps = 0;
    long int totfails = 0;
    cout << "\nRunning RKF45 with rtol = " << ARKF.rtol
         << " and atol = " << atols[ir] << endl;

    //   loop over output step sizes: call solver and output error
    while (tcur < 0.99999*Tf) {
      
      // set the time interval for this solve
      vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};

      // call the solver for this time interval
      vector<double> tvals = ARKF.Evolve(tspan, y);
      tcur = tvals.back();  // last entry in tvals
      totsteps += ARKF.steps;
      totfails += ARKF.fails;

      // compute the errors at tcur, output to screen, and accumulate maxima
      vector<double> yerr = y - ytrue(tcur);
      double abserr = InfNorm(yerr);
      double relerr = abserr / InfNorm(ytrue(tcur));
      maxrelerr = std::max(maxrelerr, relerr);
      cout << "  y(" << tcur << ") = " << setprecision(15) << y[0] 
           << endl;
      
    }

    // output final results for this tolerance
    cout << "\nOverall results for rtol = " << ARKF.rtol
     << " abstol = " << atols[ir] << ":\n"
	 << "   maxrelerr = " << maxrelerr << endl
	 << "   steps = " << totsteps << endl
	 << "   fails = " << totfails << endl;

  }
  


  // loop over relative tolerances
  for (int ir=0; ir<rtols.size(); ir++) {
	  
    // create adaptive forward Euler stepper object
    AdaptEuler AE(rhs, 0.0, atols[ir], y0);
 
    // set up the problem for this tolerance
    AE.rtol = rtols[ir];
    vector<double> y = y0;
    tcur = t0;
    double maxrelerr = 0.0;
    long int totsteps = 0;
    long int totfails = 0;
    cout << "\nRunning Adaptive Euler with rtol = " << AE.rtol
         << " and atol = " << atols[ir] << endl;

    //   loop over output step sizes: call solver and output error
    while (tcur < 0.99999*Tf) {
      
      // set the time interval for this solve
      vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};

      // call the solver for this time interval
      vector<double> tvals = AE.Evolve(tspan, y);
      tcur = tvals.back();  // last entry in tvals
      totsteps += AE.steps;
      totfails += AE.fails;

      // compute the errors at tcur, output to screen, and accumulate maxima
      vector<double> yerr = y - ytrue(tcur);
      double abserr = InfNorm(yerr);
      double relerr = abserr / InfNorm(ytrue(tcur));
      maxrelerr = std::max(maxrelerr, relerr);
      cout << "  y(" << tcur << ") = " << setprecision(15) << y[0] 
           << endl;
      
    }

    // output final results for this tolerance
    cout << "\nOverall results for rtol = " << AE.rtol 
     << " abstol = " << atols[ir] << ":\n"
	 << "   maxrelerr = " << maxrelerr << endl
	 << "   steps = " << totsteps << endl
	 << "   fails = " << totfails << endl;

  }

  return 0;
}
