/* Main routine to test Forward Euler solver for some scalar-valued ODE problems 
     y' = f(t,y), t in [0,5],
     y(0) = y0.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"
#include "fwd_euler.hpp"

using namespace std;


// ODE RHS function class -- instantiates a RHSFunction
//   includes extra routines to set up the problem and 
//   evaluate the analytical solution
class ODESystem: public RHSFunction {
  Matrix A, V, D, Vinv;
public:
  // sets up a random problem
  int Setup(int N) {
    V = Random(N,N);      // fill with random numbers
    Vinv = Inverse(V);    // Vinv = V^{-1}
    D = Matrix(N,N);
    for (int i=0; i<N; i++)
      D(i,i) = -1.2*i;
    A = V*(D*Vinv);
  }
  // evaluates the RHS function, f(t,y)
  int Evaluate(double t, vector<double>& y, vector<double>& f) {
    f = A*y;
    return 0;
  }
  // computes the true solution, y(t)
  vector<double> TrueSolution(double t, vector<double>& y0) {
    Matrix eD(y0.size(),y0.size());   // construct the matrix exponential
    for (int i=0; i<D.Rows(); i++)
      eD(i,i) = exp(D(i,i)*t);
    return (V*(eD*(Vinv*y0)));        // ytrue = V*exp(D*t)*V^{-1}*y0
  }
};



// main routine
int main(int argc, char **argv) {

  // get problem size from command line, otherwise set to 5
  int N = 5;
  if (argc > 1)
    N = atoi(argv[1]);
  cout << "\nRunning system ODE problem with N = " << N << endl;

  // time steps to try
  vector<double> h = {0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125};

  // set up problem
  ODESystem MyProblem;
  MyProblem.Setup(N);
  double t0 = 0.0;
  double Tf = 1.0;
  double dtout = 0.1;

  // initial condition
  vector<double> Y0 = Random(N);

  // create forward Euler stepper object
  ForwardEulerStepper FE(MyProblem, Y0);

  // loop over time step sizes
  for (int ih=0; ih<h.size(); ih++) {

    // set the initial condition, initial time
    vector<double> Y(Y0);
    double tcur = t0;

    // reset maxerr
    double maxabserr = 0.0;
    double maxrelerr = 0.0;

    cout << "\nRunning with stepsize h = " << h[ih] << ":\n";

    // loop over output step sizes: call solver and output error
    while (tcur < 0.99999*Tf) {
      
      // set the time interval for this solve
      vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};

      // call the solver, update current time
      vector<double> tvals = FE.Evolve(tspan, h[ih], Y);
      tcur = tvals.back();   // last entry in tvals

      // compute the error at tcur and output to screen
      vector<double> Ytrue = MyProblem.TrueSolution(tcur,Y0);
      vector<double> Yerr = Y - Ytrue;
      double abserr = InfNorm(Yerr);
      double relerr = InfNorm(Yerr) / InfNorm(Ytrue);
      cout << " Y(" << tcur << ") =" << Y;

      // accumulate overall solver statistics
      maxabserr = std::max(maxabserr, abserr);
      maxrelerr = std::max(maxrelerr, relerr);

    }
    cout << "Overall abserr = " << maxabserr 
	 << ",  relerr = " << maxrelerr << endl;
  }

  return 0;
}
