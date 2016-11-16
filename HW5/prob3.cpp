/* Main routine to test a set of LMM on the scalar-valued ODE problem 
     y' = lambda*y + 1/(1+t^2) - lambda*arctan(t), t in [0,3],
     y(0) = 0,
   using the O(h^3) AM and BDF time steppers.
   Based off of Dan Reynolds' LMM driver.
   
   Nicole Deyerl
   Math 6321 @ SMU
   Fall 2016  */

#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "rhs.hpp"
#include "trapezoidal.hpp"
#include "lmm.hpp"

using namespace std;


// Define classes to compute the ODE RHS function and its Jacobian

//    ODE RHS function class -- instantiates a RHSFunction
class MyRHS: public RHSFunction {
public:
  double lambda;                                                    // stores some local data
  int Evaluate(double t, vector<double>& y, vector<double>& f) {    // evaluates the RHS function, f(t,y)
    f[0] = lambda*y[0] + (1/(1+t*t)) - lambda*atan(t);
    return 0;
  }
};

//    ODE RHS Jacobian function class -- instantiates a RHSJacobian
class MyJac: public RHSJacobian {
public:
  double lambda;                                            // stores some local data
  int Evaluate(double t, vector<double>& y, Matrix& J) {    // evaluates the RHS Jacobian, J(t,y)
    J = 0.0;
    J(0,0) = lambda;
    return 0;
  }
};


// Convenience function for analytical solution
vector<double> ytrue(const double t) { 
  vector<double> yt = {atan(t)}; // via wolframalpha
  return yt;
};



// main routine
int main() {

  // time steps to try
  vector<double> h = {0.1, 0.01, 0.001, 0.0001};

  // lambda value
  vector<double> lambda = {-10.0, -100.0, -1000.0, -10000.0};

  // set problem information
  vector<double> y0 = {0.0};
  long int M=1;
  double t0 = 0.0;
  double Tf = 3.0;
  double dtout = 0.3;
  
  vector<double> errs(h.size());

  for (int il = 0; il<lambda.size(); il++){
	  
    // create ODE RHS and Jacobian objects, store lambda value for this test
    MyRHS rhs;
    MyJac Jac;
    rhs.lambda = lambda[il];
    Jac.lambda = lambda[il];
    
    cout << "\nlambda = " << lambda[il] << ":" << endl;

    ////////// AM-2 //////////
    // Third order Adams Moulton Method
    cout << "\n AM-2:\n";

    // create Trapezoid and LMM methods, set solver parameters
    TrapezoidalStepper TrapAM2(rhs, Jac, y0);
    vector<double> AM2_a = {1.0, 0.0};
    vector<double> AM2_b = {5.0/12.0, 8.0/12.0, -1.0/12.0};
    LMMStepper AM2(rhs, Jac, y0, AM2_a, AM2_b);
    TrapAM2.newt->SetTolerances(1.e-9, 1.e-11);
    TrapAM2.newt->SetMaxit(20);
    AM2.newt->SetTolerances(1.e-9, 1.e-11);
    AM2.newt->SetMaxit(20);

    // loop over time step sizes
    for (int ih=0; ih<h.size(); ih++) {
      cout << "  h = " << h[ih];

      // set the initial time, first output time
      double tcur = t0;
      double tout = t0 + dtout;

      // reset maxerr
      double maxerr = 0.0;
     
      // AM-2 requires two initial conditions
      Matrix y_AM2(M,2);
      //   first is just y0 (insert into 2nd column of y_AM2)
      y_AM2[1] = y0;
      //   second comes from trapezoid step step (insert into 1st column of y_AM2)
      vector<double> tspan = {tcur, tcur+h[ih]};
      y_AM2[0] = y0;
      vector<double> tvals = TrapAM2.Evolve(tspan, h[ih], y_AM2[0]);
      //   update tcur to end of initial conditions
      tcur += h[ih];

      // loop over output step sizes: call solver and output error
      while (tcur < 0.99999*Tf) {
        
        // set the time interval for this solve
        tspan = {tcur, tout};

        // call the solver, update current time
        tvals = AM2.Evolve(tspan, h[ih], y_AM2);
        tcur = tvals.back();   // last entry in tvals

        // compute the error at tcur, output to screen and accumulate maximum
        vector<double> yerr = y_AM2[0] - ytrue(tcur);   // computed solution is in 1st column of y_AM2
        double err = InfNorm(yerr);
        maxerr = std::max(maxerr, err);

        // update output time for next solve
        tout = std::min(tcur + dtout, Tf);

      }
      
      cout << "\t Max error = " << maxerr;    
      errs[ih] = maxerr;
      //convergence print out
      if (ih > 0)
        cout << "\t conv rate = " << (log(errs[ih])-log(errs[ih-1]))/(log(h[ih])-log(h[ih-1]));
          cout << endl;
      
    }
  
    ////////// BDF-3 //////////
    // Third order BDF method
    cout << "\n BDF-3:\n";

    // create Trapezoid and LMM methods, set solver parameters
    TrapezoidalStepper TrapBDF(rhs, Jac, y0);
    vector<double> BDF3_a = {18.0/11.0, -9.0/11.0, 2.0/11.0};
    vector<double> BDF3_b = {6.0/11.0, 0.0, 0.0, 0.0};
    LMMStepper BDF3(rhs, Jac, y0, BDF3_a, BDF3_b);
    TrapBDF.newt->SetTolerances(1.e-9, 1.e-11);
    TrapBDF.newt->SetMaxit(20);
    BDF3.newt->SetTolerances(1.e-9, 1.e-11);
    BDF3.newt->SetMaxit(20);

    // loop over time step sizes
    for (int ih=0; ih<h.size(); ih++) {
      cout << "  h = " << h[ih];

      // set the initial time, first output time
      double tcur = t0;
      double tout = t0 + dtout;

      // reset maxerr
      double maxerr = 0.0;
     
      // BDF-3 requires three initial conditions
      Matrix y_BDF3(M,3);
      //   first is just y0 (insert into 3rd column of y_BDF3)
      y_BDF3[2] = y0;
      //   second comes from trapezoid step step (insert into 2nd column of y_BDF3)
      vector<double> tspan = {tcur, tcur+h[ih]};
      y_BDF3[1] = y0;
      vector<double> tvals = TrapBDF.Evolve(tspan, h[ih], y_BDF3[1]);
      //   update tcur to end of initial conditions
      tcur += h[ih];
      //   third comes from trapezoid step step (insert into 1st column of y_BDF3)
      tspan = {tcur, tcur+h[ih]};
      y_BDF3[0] = y_BDF3[1];
      tvals = TrapBDF.Evolve(tspan, h[ih], y_BDF3[0]);
      //   update tcur to end of initial conditions
      tcur += h[ih];

      // loop over output step sizes: call solver and output error
      while (tcur < 0.99999*Tf) {
        
        // set the time interval for this solve
        tspan = {tcur, tout};

        // call the solver, update current time
        tvals = BDF3.Evolve(tspan, h[ih], y_BDF3);
        tcur = tvals.back();   // last entry in tvals

        // compute the error at tcur, output to screen and accumulate maximum
        vector<double> yerr = y_BDF3[0] - ytrue(tcur);   // computed solution is in 1st column of y_BDF3
        double err = InfNorm(yerr);
        maxerr = std::max(maxerr, err);

        // update output time for next solve
        tout = std::min(tcur + dtout, Tf);

      }
      
      cout << "\t Max error = " << maxerr;
      
      errs[ih] = maxerr;
      //convergence print out
      if (ih > 0)
        cout << "\t  conv rate = " << (log(errs[ih])-log(errs[ih-1]))/(log(h[ih])-log(h[ih-1]));
          cout << endl;
      
    }
  }
  return 0;
}
