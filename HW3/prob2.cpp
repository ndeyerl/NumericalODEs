/* Homework 3, problem 2: test adaptive Forward Euler for the 
   matrix ODE problem 
     y1' = -3*y1 + y2 - exp(-2*t), t in [0,3],
     y2' = y1 - 3y2 + exp(-t)
     
     matrix form:
     y' = A*y + g(t)  A=[-3,1:1,-3], g(t)=[-exp(-2*t);exp(-t)]
     
     y1(0) = 2
     y2(0) = 1

   Note: this driver uses Dr Reynolds adapt_euler code.
   Nicole Deyerl
   Math 6321 @ SMU
   Fall 2016  */

#include <iostream>
#include <iomanip>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"
#include "adapt_euler.hpp"

using namespace std;


// ODE RHS function class -- instantiates a RHSFunction
class MyRHS: public RHSFunction {
    Matrix A;
public:
  // sets up the matrix for this problem
  int Setup() {
    A = Matrix(2,2);
    A(0,0) = -3; A(0,1) = 1;
    A(1,0) = 1;  A(1,1) = -3;
  }
  // sets up the vector g(t) and the rhs f
  int Evaluate(double t, vector<double>& y, vector<double>& f) {
    vector<double> g(2);
    g[0] = -exp(-2.0*t);
    g[1] = exp(-t);
    f = A*y + g;
    return 0;
  }
};
// Analytic solution to the ode system
vector<double> ytrue(double t) { 
    vector<double> eD(2);   
    eD[0] = (1.0/3.0)*exp(-t) + ((3.0/4.0)-t/2.0)*exp(-2.0*t)  
                + (11.0/12.0)*exp(-4.0*t);
    eD[1] =  (2.0/3.0)*exp(-t) + ((5.0/4.0)-t/2.0)*exp(-2.0*t)  
                - (11.0/12.0)*exp(-4.0*t);
    return eD;   
};


// main routine
int main() {

  // tolerances to try
  vector<double> rtols = {1.e-2, 1.e-4, 1.e-6, 1.e-8};
  double atol = 1.e-12;
  
  // initial condition and time span given by problem
  vector<double> y0 = {2.0,1.0}; // vector containing initial condition
  double t0 = 0.0;
  double Tf = 3.0;
  double tcur = t0;
  double dtout = 0.3;

  // create ODE RHS function objects
  MyRHS rhs;
  rhs.Setup();

  // create forward Euler stepper object (will reset rtol before each solve)
  AdaptEuler AE(rhs, 0.0, atol, y0);

  // loop over relative tolerances
  for (int ir=0; ir<rtols.size(); ir++) {
 
    // set up the problem for this tolerance
    AE.rtol = rtols[ir];
    vector<double> y = y0;
    tcur = t0;
    double maxabserr = 0.0; //initialize holder for the error between dtouts
    double maxy = 0.0; // initialize holder over all y-values (for rtol*||y||...)
    long int totsteps = 0;
    long int totfails = 0;
    cout << "\nRunning problem 2 with rtol = " << AE.rtol
         << " and atol = " << atol << endl;

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
      maxy = std::max(maxy, InfNorm(y));
      maxabserr = std::max(maxabserr, abserr);
      cout << "  t = " << tcur << "  y1 = " << setprecision(17) << y[0] 
           << "  y2 = " << setprecision(17) << y[1]
           << setprecision(6)
           << "\t abserr = " << abserr 
           << endl;
      
    }

    // output final results for this tolerance
    cout << "\nOverall results for rtol = " << AE.rtol << ":\n"
         << "   max-norm err = " << maxabserr << endl
         << "   rtol*||y||+atol = " << rtols[ir]*maxy + atol << endl
	     << "   steps = " << totsteps << endl;

  }

  return 0;
}
