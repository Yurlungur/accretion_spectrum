// bisection.hpp

// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2013-12-07 22:15:29 (jonah)>

// This is prototype for a general bisection root-finding algorithm
// for use in shooting method applications.

// Include guard
# pragma once
// Includes
#include <vector>   // for output and internal variables
#include "rkf45.hpp" // For the Runge-Kutta integrator
#include <float.h> // for machine epsilon
#include <cmath> // for math


// The root-finder error.
const double ROOT_ERROR = 1E10 * DBL_EPSILON;
// Are we debugging
const bool DEBUGGING = false;

// ScalarFunction class
// ----------------------------------------------------------------------

// We use this class to store the function that the root finder finds
// the root of. When called as a "function" (using the overloaded ()
// operator), it takes a double and returns a double. Other than that,
// implementation details are left to the subclass.

// Base class
class ScalarFunction {
public:
  // The "input"/"output" of the function
  virtual double operator()(double) =0;
  // times_called keeps track of how many times a functor is called.
  int get_times_called() const { return times_called; }

protected:
  int times_called;
  ScalarFunction() { times_called = 0; }
};

// We define two subclasses. One is just for a basic function. This is
// for unit testing. The other subclass is for the shooting method. It
// holds an RKF45 class internally.

// Basic function
class BasicFunction: public ScalarFunction {
public:
  // Constructor. Takes a function y(x) as input. 
  BasicFunction(double (*y)(double));
  // Destructor
  ~BasicFunction() {}

  // The input/output
  double operator()(double x);
protected:
  // The function
  double (*y)(double);
};
// ----------------------------------------------------------------------


// A basic binomial to test root finding
double binomial(double x);

// The root finder method.
// ----------------------------------------------------------------------

// Takes a functor (a bound on x) which represents a function to find
// the root of. Also takes a minimum bound on the variable x and a
// maximum bound. Returns the value for x which is the zero fo the
// ScalarFunction.The functor is explicitly allowed to be modified, as at the
// end of the day, we want the RKF45 integrator to have integrated to
// the right spot.
double bisection_root_finder(ScalarFunction& f, double x_min, double x_max);

// ----------------------------------------------------------------------
