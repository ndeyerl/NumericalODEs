/* Homework 3, problem 5: test generalized trapezoid for the 
   scalar ODE problem 
     y' = -lambda*y + 1/(1+t^2) - lambda*atan(t), t in [0,1]
     y(0) = 0
   for various stiffness parameters lambda, step sizes h and 
   trapezoid method parameters theta.
   
   Note: this driver was based off of Dr Reynolds Bwd_euler driver 
   and trapezoid cpp and hpp files, and uses Dr Reynolds Newton  
   and residual files.
   Nicole Deyerl
   Math 6321 @ SMU
   Fall 2016  */

#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"
#include "trapezoidal.hpp"

using namespace std;


// Define classes to compute the ODE RHS function and its Jacobian

//    ODE RHS function class -- instantiates a RHSFunction
class MyRHS: public RHSFunction {
public:
  double lambda;                                    // stores some local data
  int Evaluate(double t, vector<double>& y, vector<double>& f) {    // evaluates the RHS function, f(t,y)
    f[0] = lambda*y[0] + (1/(1+pow(t,2))) - lambda*atan(t); // given by problem
    return 0;
  }
};

//    ODE RHS Jacobian function class -- instantiates a RHSJacobian
class MyJac: public RHSJacobian {
public:
  double lambda;                                            // stores some local data
  int Evaluate(double t, vector<double>& y, Matrix& J) {    // evaluates the RHS Jacobian, J(t,y)
    J(0,0) = lambda; // from differentiating f given by problem
    return 0;
  }
}; 


// Convenience function for analytical solution
vector<double> ytrue(const double t) { 
  vector<double> yt = {atan(t)}; // given by problem
  return yt;
};



// main routine
int main() {

  // time steps to try
  vector<double> h = {0.1, 0.01, 0.001, 0.0001};

  // lambda values to try
  vector<double> lambdas = {-200.0, -2000.0, -20000.0};
  
  // theta values to try
  vector<double> thetas = {1.0, 0.55, 0.5, 0.45};
  
  // set problem information
  vector<double> y0 = {0.0}; // initial condition
  double t0 = 0.0;
  double Tf = 1.0;
  double dtout = 0.1;

  // create ODE RHS and Jacobian objects
  MyRHS rhs;
  MyJac Jac;
 

  //------ Trapezoidal tests ------

  // loop over theta values
  for (int it = 0; it<thetas.size(); it++) {

    // create time stepper objects
    // theta is a part of the object, so need to make a new one after each
    //loop over theta
    TrapezoidalStepper Tr(rhs, Jac, thetas[it], y0);

    // update Newton solver parameters
    Tr.newt->SetTolerances(1.e-11, 1.e-13);
    Tr.newt->SetMaxit(50);
  
    // loop over lambda values
    for (int il=0; il<lambdas.size(); il++) {
    
      // set current lambda value into rhs and Jac objects
      rhs.lambda = lambdas[il];
      Jac.lambda = lambdas[il];
    
      //error storage
      vector<double> abserrs(h.size());
    
      // loop over time step sizes
      for (int ih=0; ih<h.size(); ih++) {
      
        // set the initial condition, initial time
        vector<double> y(y0);
        double tcur = t0;

        // reset maxerr
        double maxerr = 0.0;
     
        cout << "\nRunning trapezoidal with stepsize h = " << h[ih] 
        << ", lambda = " << lambdas[il] << ", theta = " << thetas[it] << ":\n";

        // loop over output step sizes: call solver and output error
        while (tcur < 0.99999*Tf) {
      
          // set the time interval for this solve
          vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};

          // call the solver, update current time
          vector<double> tvals = Tr.Evolve(tspan, h[ih], y);
          tcur = tvals.back();   // last entry in tvals

          // compute the error at tcur, output to screen and accumulate maximum
          vector<double> yerr = y - ytrue(tcur);
          double err = InfNorm(yerr);
          maxerr = std::max(maxerr, err);
          cout << "  y(" << tcur << ") = " << y[0]
          << "  \t||error|| = " << err
          << endl;
        }
        abserrs[ih] = maxerr;
        cout << "Max error = " << maxerr << endl;
      
      }
    
      // calculate orders of convergence between successive values of h (absolute error)
      cout << "\nConvergence order estimates:\n";
      for (int ih=0; ih<h.size()-1; ih++) {
        double dlogh = log(h[ih+1]) - log(h[ih]);
        double dloge = log(abserrs[ih+1]) - log(abserrs[ih]);
        cout << "  h = " << h[ih] << "  order = " << dloge/dlogh << endl;
      }
    }
  }
  return 0;
}
