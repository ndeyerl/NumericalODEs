
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

  // time steps to try
  vector<double> h = {0.2,0.1,0.05,0.025,0.0125}; //specified by 
                                                  // problem 1a

  // set problem information
  vector<double> y0_1 = {1.0}; //initial condition y(0) = 1
  double t0 = 0.0;
  double Tf = 5.0; //problem specified t in [0,5]
  double tcur = t0;
  double dtout = 1.0;//problem specified output of soln and abs err every
                     // 1 unit of time
  double convg = 0.0; //initialize order of convergence
  double diff;
  // create ODE RHS function objects
  RHS1 f1;
  vector<double> maxerr1 (h.size(),0.0); //preallocate space to save all maxerrors
                                        // for each step size h, populate w/ 0's
  
  // create forward Euler stepper object
  ForwardEulerStepper FE1(f1, y0_1);
 

  // loop over time step sizes
  for (int ih=0; ih<h.size(); ih++) {

    // problem 1:
    vector<double> y = y0_1;
    tcur = t0;
    double maxerr = 0.0;

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
      double err = InfNorm(yerr); //absolute error = diff between num + true soln
                                  // (InfNorm -> gives maximal entry)
      maxerr1[ih] = err; //keep the errors to calculate convergence
      maxerr = std::max(maxerr, err); //keep maximal error value
      // soln + error print out
      cout << "  y(" << tcur << ") = " << y[0]
	   << "  \t||error|| = " << err
	   << endl;
    }
    cout << "Max error = " << maxerr << endl;
    if(ih > 0){
      convg = (log(maxerr1[ih])-log(maxerr1[ih-1]))/(log(h[ih])-log(h[ih-1])); //soln of err = h^p for p
      cout << "The order of convergence is = " << convg << endl;
    }
  }
  

  return 0;
}
