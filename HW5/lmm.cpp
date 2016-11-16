/* Linear multistep time stepper class implementation file.

   Class to perform time evolution of the IVP
        y' = f(t,y),  t in [t0, Tf],  y(t0) = y0
   using a linear multistep time stepping method.  Although this 
   class is written to directly support implicit LMM, it will work
   equally well for explicit LMM.  Also, the initial condition fed 
   in to Newton solver is found using a step of explicit forward Euler.
   Based off of Dan Reynolds lmm.cpp, which was just slightly altered.
   
   Nicole Deyerl
   Math 6321 @ SMU
   Fall 2016  */

#include "lmm.hpp"



//////////// LMM Residual ////////////

// residual initialization routine
int LMMResid::Initialize(double t, double h, Matrix &y) {

  // copy y into stored yold object
  int ierr = yold.Copy(y);
  if (ierr != 0) {
    std::cerr << "Error in LMMResid::Initialize Copy call = " << ierr << "\n";
    return ierr;
  }
  
  // fill initial set of 'old' RHS vectors
  for (int j=0; j<y.Columns(); j++) {
    ierr = frhs->Evaluate(t-j*h, yold[j], fold[j]);
    if (ierr != 0) {
      std::cerr << "Error in LMMResid::Initialzie ODE RHS Evaluate call = " << ierr << "\n";
      return ierr;
    }
  }

  // return success
  return 0;
}

// residual evaluation routine
int LMMResid::Evaluate(std::vector<double>& y, std::vector<double>& resid) {

  // evaluate RHS function at new time (store in resid)
  int ierr = frhs->Evaluate(t+h, y, resid);
  if (ierr != 0) {
    std::cerr << "Error in ODE RHS function = " << ierr << "\n";
    return ierr;
  }

  // combine pieces to fill residual, y - sum[a_j y_{n-j}] - h*sum[b_j*f_{n-j}]
  resid *= (-h*b[0]);              // resid = -h*b_{-1}*f(t+h,y)
  resid += y;                      // resid = y - h*b_{-1}*f(t+h,y)
  for (int j=0; j<a.size(); j++)
    resid -= ( a[j] * yold[j] + (h*b[j+1]) * fold[j] );
  
  // return success
  return 0;
}

// Routine to handle updates of "old" solutions and right-hand sides
//
// Inputs:  tnew  the current time for the new solution
//          ynew  the new solution
// Returns: 0 (success) or 1 (failure)
int LMMResid::Update(double t, std::vector<double>& ynew) {

  // update columns of yold and fold, starting at oldest and moving to newest
  for (int icol=fold.Columns()-1; icol>0; icol--) {
    fold[icol] = fold[icol-1];
    yold[icol] = yold[icol-1];
  }

  // fill first column of yold with ynew
  yold[0] = ynew;

  // evaluate RHS function at new time (store in first column of fold matrix)
  int ierr = frhs->Evaluate(t, ynew, fold[0]);
  if (ierr != 0) {
    std::cerr << "Error in ODE RHS function = " << ierr << "\n";
    return ierr;
  }

  return 0;
};


//////////// LMM Residual Jacobian ////////////

// Jacobian evaluation routine
int LMMResidJac::Evaluate(std::vector<double>& y, Matrix& J) {

  // evaluate RHS function Jacobian (store in J)
  int ierr = Jrhs->Evaluate(t+h, y, J);
  if (ierr != 0) {
    std::cerr << "Error in ODE RHS Jacobian function = " << ierr << "\n";
    return ierr;
  }
  // combine pieces to fill residual Jacobian
  J *= (-beta*h);                   // J = -beta*h*Jrhs
  for (int i=0; i<J.Rows(); i++)    // J = I - beta*h*Jrhs
    J(i,i) += 1.0;

  // return success
  return 0;
}


//////////// LMM Time Stepper ////////////

// LMM construction routine (allocates local data)
//
// Inputs:  frhs_  holds the RHSFunction to use
//          Jrhs_  holds the RHSJacobian to use
//          y      holds a template solution vector (for cloning)
//          a,b    hold the LMM coefficients
LMMStepper::LMMStepper(RHSFunction& frhs_, 
                       RHSJacobian& Jrhs_, 
		       std::vector<double>& y, 
                       std::vector<double>& a_, 
                       std::vector<double>& b_) {

  // check that LMM coefficients are compatible
  if (b_.size() != a_.size()+1) {
    std::cerr << "LMMStepper: Incompatible LMM coefficients; will not function!\n";
    return;
  }

  // point at and/or copy inputs
  frhs = &frhs_;   // set pointer to RHSFunction
  a = a_;          // copy LMM coefficient arrays
  b = b_;

  // construct objects for nonlinear residual and its Jacobian
  resid = new LMMResid(frhs_, y, a, b);
  residJac = new LMMResidJac(Jrhs_, b[0]);

  // construct Newton solver object (only copies y)
  // (initialize with default solver parameters; user may override)
  newt = new NewtonSolver(*resid, *residJac, 1.0e-7, 1.0e-11, 100, y, false);
};


// The actual LMM time step evolution routine
//
// Inputs:  tspan holds the current time interval, [t0, tf]
//          h holds the desired time step size
//          y holds the set of initial conditions, [y(0), y(-1), ..., y(-p)]
// Outputs: y holds the computed solution and set of p previous solutions, 
//             [y(tf), y(tf-h), ..., y(tf-p*h)]
//
// The return value is a row vector containing all internal 
// times at which the solution was computed,
//               [t0, t1, ..., tN]
std::vector<double> LMMStepper::Evolve(std::vector<double>& tspan, 
                                       double h, Matrix& y) {

  // initialize output
  std::vector<double> times = {tspan[0]};

  // check for legal inputs 
  if (h <= 0.0) {
    std::cerr << "Evolve: Illegal h\n";
    return times;
  }
  if (tspan[1] <= tspan[0]) {
    std::cerr << "Evolve: Illegal tspan\n";
    return times;	  
  }
  if (y.Columns() != a.size()) {
    std::cerr << "LMMStepper: insufficient set of initial conditions, y.Columns() = " 
	      << y.Columns() << " and a.size() = " << a.size() << "\n";
    return times;
  }
  
  // initialize residual structures 
  int ierr = resid->Initialize(tspan[0], h, y);
  if (ierr != 0) {
    std::cerr << "LMMStepper::Evolve error in residual Initialize call = " << ierr << "\n";
    return times;
  }
  
  std::vector<double> f = y[0]; //create a vector to hold the feval
  
  // compute ODE RHS
  if (frhs->Evaluate(times[0], y[0], f) != 0) {
    std::cerr << "ForwardEulerStepper: Error in ODE RHS function\n";
    return times;
  }
  
  // set vector ycur to contain one step of explicit euler
  std::vector<double> ycur = y[0];

  // figure out how many time steps
  long int N = (tspan[1]-tspan[0])/h;
  if (tspan[1] > tspan[0]+N*h)  N++;

  // iterate over time steps
  for (int i=0; i<N; i++) {

    // last step only: update h to stop directly at final time
    if (i == N-1) 
      h = tspan[1]-times[i];

    // update resid and residJac objects with information on current state
    resid->t    = times[i];     // copy current time into objects
    residJac->t = times[i];
    resid->h    = h;            // copy current stepsize into objects
    residJac->h = h;

    ycur += h*f; // one step of forward euler per time step -> one IC per time step
    
    // call Newton method to solve for the updated solution
    // (initial guess is one step of forward euler)
    int nerr = newt->Solve(ycur);
    if (nerr != 0) {
      std::cerr << "LMMStepper: Error in Newton solver function = " << nerr << "\n";
      return times;
    }

    // update current time (store in output array), and "old" data
    times.push_back(times[i] + h);
    resid->Update(times[i+1], ycur);
  }

  // copy yold data back into y Matrix
  y.Copy(resid->yold);

  return times;
}
