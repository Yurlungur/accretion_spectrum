// bisection.cpp

// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2013-12-04 16:42:55 (jonah)>

// This is the implementation for a general bisection root-finding
// algorithm for use in shooting method applications.

// Icludes
#include "bisection.hpp"
#include <iostream> // for debugging
#include <cassert>
using std::abs;
using std::cout;
using std::endl;

// Functors
// ------------------------------------------------------------
double BasicFunction::operator()(double x) {
  times_called++;
  return y(x);
}
BasicFunction::BasicFunction(double (*y)(double)) {
  this->y = y;
}
// ----------------------------------------------------------------------


// Utilities
// ----------------------------------------------------------------------
double binomial(double x) {
  const double Y_OFFSET = 16;
  return x*x - Y_OFFSET;
}
// ----------------------------------------------------------------------


// The root finder method
// ----------------------------------------------------------------------
double bisection_root_finder(ScalarFunction& f, double x_min, double x_max) {
  assert ( x_min < x_max );
  if (DEBUGGING) {
    cout << "x_min = " << x_min << endl;
    cout << "x_max = " << x_max << endl;
  }
  double rightmost_y = f(x_max);
  if (DEBUGGING) {
    cout << "\trightmost_y = " << rightmost_y << endl;
  }
  if ( abs(rightmost_y) < ROOT_ERROR ) { // Then we're done!
    return x_max;
  }
  bool f_increasing = rightmost_y > 0;
  double midpoint = x_min + 0.5 * (x_max - x_min);
  double mid_y = f(midpoint);
  if (DEBUGGING) {
    cout << "\tmidpoint = " << midpoint << endl;
    cout << "\tmid_y = " << mid_y << endl;
  }
  if ( abs(mid_y) < ROOT_ERROR ) { // then we're done!
    return midpoint;
  }
  else if ( ( f_increasing && mid_y > 0 ) || (!(f_increasing) && mid_y < 0) ) {
    // Then we're too far right
    return bisection_root_finder(f,x_min,midpoint);
  }
  else {
    // We're too far left
    return bisection_root_finder(f,midpoint,x_max);
  }
}
