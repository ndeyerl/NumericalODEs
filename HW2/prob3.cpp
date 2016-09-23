
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
#include "adapt_euler.hpp"

using namespace std;


// Problem 1:  y' = -e^(-t)*y, y(0)=1
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



// main routine
int main() {

  // tolerances to try
  vector<double> rtol = {pow(10,-2),pow(10,-4),pow(10,-6),pow(10,-8)};
  double atol = pow(10,-11);

  // set problem information
  vector<double> y0_1 = {1.0}; //initial condition y(0) = 1
  double t0 = 0.0;
  double Tf = 5.0; //problem specified t in [0,5]
  double tcur = t0;
  double dtout = 1.0;//problem specified output of soln and abs err every
                     // 1 unit of time
  double diff;
  // create ODE RHS function objects
  RHS1 f1;
  vector<double> maxrel (rtol.size(),0.0); //preallocate space to save all maxerrors
                                        // for each step size h, populate w/ 0's
  vector<double> maxabs (rtol.size(),0.0);

  // loop over time step sizes
  for (int ir=0; ir<rtol.size(); ir++) {

    // create an adaptive Euler stepper object for each rtol
    AdaptEuler FE1(f1, rtol[ir], atol, y0_1);
  
    // problem 1:
    vector<double> y = y0_1;
    tcur = t0;
    double maxabserr = 0.0;
    double maxrelerr = 0.0;
    double relerr = 0.0;
    double abserr = 0.0;

    cout << "\nRunning problem 1 with rtol = " << rtol[ir] << " atol = " << atol << ":\n";

    //   loop over output step sizes: call solver and output error
    while (tcur < 0.9999*Tf) {
      
      // set the time interval for this solve
      vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};
      
      // call the solver, update current time
      vector<double> tvals = FE1.Evolve(tspan, y);
      tcur = tvals.back();   // last entry in tvals
      // compute the error at tcur, output to screen and accumulate maximum
      vector<double> yerr = y - ytrue1(tcur);
      double abserr = InfNorm(yerr);
      double relerr = InfNorm(yerr)/InfNorm(y);
      maxabserr = std::max(maxabserr, abserr);
      maxrelerr = std::max(maxrelerr, relerr);
      cout << "  y(" << tcur << ") = " << y[0]
	   << "  \t||abs. error|| = " << abserr << "  \t||rel. error|| = " << relerr
	   //<< "  \t||error|| = " << err
	   << endl;
    }
    cout << "Number of calls to f = " << FE1.fcalls << std::endl;
    cout << "Max absolute error = " << maxabserr 
         << "  Max relative error = " << maxrelerr << endl;

  }
  

  return 0;
}
