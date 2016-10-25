/* Testing routine for generalized trapezoid solver from homework 3

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016  */

#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"
#include "trapezoidal.hpp"

using namespace std;


// Define classes to compute the ODE RHS function and its Jacobian
//   ODE RHS function class -- instantiates a RHSFunction
//   includes extra routines to set up the problem and 
//   evaluate the analytical solution
class MyRhs: public RHSFunction {
public:
  Matrix A, V, D, Vinv;
  // sets up a random problem
  int Setup(int N) {
    V = Random(N,N);      // fill with random numbers
    Vinv = Inverse(V);    // Vinv = V^{-1}
    D = Matrix(N,N);
    for (int i=0; i<N; i++)
      D(i,i) = -(500.0*(i+1))/N;
    A = V*(D*Vinv);
    return 0;
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

//    ODE RHS Jacobian function class -- instantiates a RHSJacobian
class MyJac: public RHSJacobian {
public:
  Matrix A;                      // stores some local data
  int Evaluate(double t, vector<double>& y, Matrix& J) {    // evaluates the RHS Jacobian, J(t,y)
    J = A;
    return 0;
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
  vector<double> h = {0.01, 0.005, 0.002, 0.001, 0.0005};

  // theta values to try
  vector<double> thetas = {0.55, 0.5, 0.45, 0.35, 0.25, 0.15};

  // set up problem, including ODE RHS and Jacobian objects
  vector<double> y0 = 50.0*Random(N);
  double t0 = 0.0;
  double Tf = 1.0;
  double dtout = 0.1;
  MyRhs rhs;
  rhs.Setup(N);
  MyJac Jac;
  Jac.A = rhs.A;

  cout << "Running Generalized Trapezoid tests:\n";


  // loop over theta values
  for (int it=0; it<thetas.size(); it++) {

    // create time stepper object; update Newton solver parameters
    //GeneralizedTrapezoid GT(rhs, Jac, thetas[it], y0);
    TrapezoidalStepper GT(rhs, Jac, thetas[it], y0);
    GT.newt->SetTolerances(1.e-14, 1.e-15);
    GT.newt->SetMaxit(20);

    // create array for error values
    vector<double> errors(h.size());

    // loop over time step sizes
    for (int ih=0; ih<h.size(); ih++) {
      
      // set the initial condition, initial time
      vector<double> y(y0);
      double tcur = t0;

      // reset maxerr
      double maxerr = 0.0;
     
      cout << "  theta = " << thetas[it]
	   << ",  \th = " << h[ih] << ":";

      // loop over output step sizes: call solver and output error
      while (tcur < 0.99999*Tf) {
      
	// set the time interval for this solve
	vector<double> tspan = {tcur, std::min(tcur + dtout, Tf)};

	// call the solver, update current time
	vector<double> tvals = GT.Evolve(tspan, h[ih], y);
	tcur = tvals.back();   // last entry in tvals

	// compute the error at tcur, output to screen and accumulate maximum
	vector<double> ytrue = rhs.TrueSolution(tcur,y0);
	vector<double> yerr = y - ytrue;
	double err = InfNorm(yerr);
	maxerr = std::max(maxerr, err);
      }
      cout << " \terror = " << maxerr << endl;
      errors[ih] = maxerr;

    }

    // output convergence rate estimates
    cout << "  Convergence rates:\n  ";
    for (int ih=1; ih<h.size(); ih++)
      cout << "  \t" << (log(errors[ih])-log(errors[ih-1]))/(log(h[ih])-log(h[ih-1]));
    cout << endl << endl;

  }

  return 0;
}
