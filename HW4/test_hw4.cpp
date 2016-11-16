/* Homework 4 testing program.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"
#include "ab3.hpp"
#include "erk4.hpp"
#include "adapt_rkf.hpp"

using namespace std;


// ODE RHS function class -- instantiates a RHSFunction
//   includes extra routine to evaluate the analytical solution
class ODESystem: public RHSFunction {
public:
  double lambda;
  // evaluates the RHS function, f(t,y)
  int Evaluate(double t, vector<double>& y, vector<double>& f) {
    f[0] = (6.0+4.0*lambda)*y[0] - (4.0+2.0*lambda)*y[1] - (3.0+2.0*lambda)*y[2] + 5.0*exp(-t);
    f[1] = (6.0*lambda-3.0)*y[0] + (2.0-3.0*lambda)*y[1] + (2.0-3.0*lambda)*y[2] + 6.0*exp(-t);
    f[2] = 15.0*y[0] - 10.0*y[1] - 8.0*y[2] + 2.0*exp(-t);
    return 0;
  }
  // computes the true solution, y(t)
  vector<double> TrueSolution(double t) {
    vector<double> yt(3);
    yt[0] = -(9.0+lambda)/(2.0*lambda+2.0)*exp(-t)
          - 4.0*lambda/(lambda+1.0)*exp(lambda*t)
          + 2.5*cos(t) + 0.5*sin(t);
    yt[1] = (lambda-11.0)/(2.0*lambda+2.0)*exp(-t)
          - 6.0*lambda/(lambda+1.0)*exp(lambda*t) 
          - 0.5*cos(t) + 2.5*sin(t);
    yt[2] = -1.5*exp(-t) + 5.5*cos(t) - 1.5*sin(t);
    return yt;
  }
};



// main routine
int main() {

  // set up problem
  vector<double> y0 = {-2.0, -6.0, 4.0};
  double t0 = 0.0;
  double Tf = 5.0;
  double tcur = t0;
  double dtout = 0.2;
  double lambda = -20.0;

  // create ODE RHS function object
  ODESystem f;
  f.lambda = lambda;
 

  //----- first test ab3 -----//

  // time steps to try
  vector<double> h = {0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002};

  // create forward Euler stepper object
  AB3Stepper AB3(f, y0);

  // create RK4 stepper object
  ERK4Stepper RK4(f, y0);

  // storage for errors
  vector<double> abserrs(h.size());

  cout << "\nAB3 time stepper:\n";

  // loop over time step sizes
  for (int ih=0; ih<h.size(); ih++) {

    cout << "  h = " << h[ih];

    // initialize maximum error 
    double tout = t0+dtout;
    double maxerr = 0.0;
    
    // generate extra initial conditions via ERK4
    tcur = t0;
    vector<double> y3 = y0;
    vector<double> y2 = y3;
    RK4.Step(tcur, h[ih], y2);
    tcur += h[ih];
    vector<double> y1 = y2;
    RK4.Step(tcur, h[ih], y1);
    tcur += h[ih];
    vector<double> y = y1;
    RK4.Step(tcur, h[ih], y);
    tcur += h[ih];

    //   loop over output step sizes: call solver and output error
    while (tcur < 0.99999*Tf) {
      
      // set the time interval for this solve
      vector<double> tspan = {tcur, std::min(tout, Tf)};

      // call the solver for this time interval
      vector<double> tvals = AB3.Evolve(tspan, h[ih], y, y1, y2, y3);
      tcur = tvals.back();   // last entry in tvals
      tout = tcur + dtout;

      // compute the error at tcur, output to screen and accumulate maximum
      vector<double> yerr = y - f.TrueSolution(tcur);
      maxerr = std::max(maxerr, InfNorm(yerr));
    }
    abserrs[ih] = maxerr;
    cout << ",\t error = " << abserrs[ih] << "\n";

  }

  // calculate orders of convergence between successive values of h (absolute error)
  cout << "\nAB3 convergence order estimates:\n";
  for (int ih=0; ih<h.size()-1; ih++) {
    double dlogh = log(h[ih+1]) - log(h[ih]);
    double dloge = log(abserrs[ih+1]) - log(abserrs[ih]);
    cout << "  order = " << dloge/dlogh << endl;
  }


  //----- second test adapt_rkf -----//

  // tolerances
  vector<double> rtols = {1.e-6, 1.e-8, 1.e-10};
  vector<double> atols = {1.e-8, 1.e-10, 1.e-12};

  // storage for errors at each tolerance
  vector<double> err_RKF(rtols.size());


  /////////////  Runge-Kutta-Fehlberg /////////////
  cout << "\nAdaptive RKF solver, steps and errors vs tolerances:\n";

  // loop over tolerances
  for (int ir=0; ir<rtols.size(); ir++) {

    // create adaptive RKF solver (will reset tolerances)
    AdaptRKF RKF(f, rtols[ir], atols[ir], y0);

    // set the initial condition, initial times
    vector<double> y = y0;
    double tcur = t0;

    // reset maxerr
    double maxerr = 0.0;

    cout << "  rtol = " << rtols[ir] << ",  atol = " << atols[ir];
	
    // loop over output steps: call solver and output error
    while (tcur < 0.99999*Tf) {

      // set the time interval for this call
      std::vector<double> tspan = { tcur, std::min( tcur + dtout, Tf) };

      // call the solver, update current time
      vector<double> tvals = RKF.Evolve(tspan, y);
      tcur = tvals.back();   // last entry in tvals

      // compute the error at tcur and accumulate statistics
      vector<double> yerr = y - f.TrueSolution(tcur);
      maxerr = std::max(maxerr, InfNorm(yerr));
    }
    cout << ",  error = " << maxerr << ",  steps = " 
	 << RKF.steps << ",  fails = " << RKF.fails << endl;
  }

  return 0;
}
