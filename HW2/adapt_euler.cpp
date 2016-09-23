/* Forward Euler time stepper class implementation file.

   Class to perform time evolution of the IVP
        y' = f(t,y),  t in [t0, Tf],  y(t0) = y0
   using the forward Euler (explicit Euler) time stepping method. 

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2016 */

#include <vector>
#include "matrix.hpp"
#include "adapt_euler.hpp"


// The forward Euler time step evolution routine
//
// Inputs:  tspan holds the current time interval, [t0, tf]
//          h holds the desired time step size
//          y holds the initial condition, y(t0)
// Outputs: y holds the computed solution, y(tf)
//
// The return value is a row vector containing all internal 
// times at which the solution was computed,
//               [t0, t1, ..., tN]
std::vector<double> AdaptEuler::Evolve(std::vector<double> tspan, std::vector<double>& y) {

  // initialize output
  std::vector<double> times = {tspan[0]};

  // check for legal inputs 
  if (tspan[1] <= tspan[0]) {
    std::cerr << "AdaptEuler: Illegal tspan\n";
    return times;	  
  }

  // figure out how many time steps
  long int N = pow(10,6);
  double h = (tspan[1]-tspan[0])/(tspan[1]*8);
//  double h = .001;
  double err = 0.0;
  double relerr, abserr;
  double th, th2;
  int counter = 0; //counter for index of time steps
  std::vector<double> yh = y;
  std::vector<double> yh2 = y;

  // iterate over time steps
  while (times[counter]<0.99999*tspan[1]) {
	
	//set up times for evals (at half and full step) 
    th = times[counter] + h;
    th2 = times[counter] + h/2.0;
    
    //eval the full step
    if (frhs->Evaluate(th, y, f) != 0) {
      std::cerr << "AdaptEuler: Error in ODE RHS function\n";
      return times;
    }
    
    //1 full euler step, 1 half euler step
    yh = y + h*f; //1 full step
    yh2 = y + (h/2.0)*f; //2 half steps
    
    //eval the half step
    if (frhs->Evaluate(th2, yh2, f) != 0) {
      std::cerr << "AdaptEuler: Error in ODE RHS function\n";
      return times;
    }
    
    //second half euler step
    yh2 += (h/2.0)*f;
    
    //richardson error estimate
    err = Norm(2.0*yh2 - 2.0*yh);
    
    
    if(err <= r*Norm(y) + a){
		y = 2.0*yh2 - yh; //richardson euler formula
		counter = counter + 1; //update counter
        times.push_back(th); //update current time, store in output array
        fcalls += 2; //update number of calls to f
	} else {
		h = h*((r*Norm(y) + a)/err); //reduce step size + try again
		fcalls += 2; //update number of calls to f
	}
    
    if(counter >= pow(10,6)){ //if the number of time steps exceeds 10^6
		break;                // exit the solver
	}

  }
  return times;
}

