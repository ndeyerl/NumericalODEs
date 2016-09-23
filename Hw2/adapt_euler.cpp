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
  double h = .01;
  double err = 0.0;
  int counter = 0;
  std::vector<double> yfull = y;
  std::vector<double> yhalf = y;

  // iterate over time steps
  for (long int i=0; i<N; i++) {

    while(times[counter]<tspan[1]){
    // compute ODE RHS
    if (frhs->Evaluate(times[counter], y, f) != 0) {
      std::cerr << "AdaptEuler: Error in ODE RHS function\n";
      return times;
    }
    if (h*counter > tspan[1]){
		break;
	}

    // update solution with forward Euler step

    yfull = y + (h*f); //1 full step
    yhalf = y + ((h/2.0)*f); //2 half steps
    //yhalf += (h*f);

    if(err <= r*2.0*yhalf[counter] + a){
		
		y = 2.0*yhalf - yfull;
		h = h*err;
        times.push_back((counter+1)*h); // update current time, store in output array
        std::cout << "time = " << y << std::endl;
        counter = counter + 1;
	}else{
		h = h/err;

	}
    

  }
  std::cout << "times" << std::endl;
  return times;

}
}
