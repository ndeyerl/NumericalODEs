/* Explicit 4th order Adams-Bashforth time stepper class header file.
 * Based on Dan Reynolds ab1.hpp file.

   Nicole Deyerl
   Math 6321 @ SMU
   Fall 2016  */

#ifndef AB3_DEFINED__
#define AB3_DEFINED__

// Inclusions
#include <math.h>
#include <vector>
#include "matrix.hpp" 
#include "rhs.hpp"


// Explicit AB3 time stepper class
class AB3Stepper {

 private:

  RHSFunction *frhs;            // pointer to ODE RHS function
  std::vector<double> f, f1, f2, f3;  // reused vectors
  
 public:

  // constructor (sets RHS function pointer, allocates local data)
  AB3Stepper(RHSFunction& frhs_, std::vector<double>& y) { 
    frhs = &frhs_;      // store RHSFunction pointer
    f = y;              // allocate reusable data
    f1 = y;             // based on size of y, holds old f (f at yn-1)
    f2 = y;             // holds older f (f at yn-2)
    f3 = y;             // holds oldest f (f at yn-3)
  };

  // Evolve routine (evolves the solution)
  std::vector<double> Evolve(std::vector<double>& tspan, 
			     double h, 
                 std::vector<double>& y, 
			     std::vector<double>& y1,
			     std::vector<double>& y2,
			     std::vector<double>& y3);

};

#endif
