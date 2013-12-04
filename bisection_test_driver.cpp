// bisection_test_driver.cpp

// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2013-12-04 16:37:25 (jonah)>

// This tests the bisection root finding algorithm.

// Includes
#include "bisection.hpp"
#include <iostream>
using std::cout;
using std::endl;

const double XMIN = 0;
const double XMAX = 10;

int main() {
  cout << "Testing the bisection root finder algorithm." << endl;
  cout << "Building functor." << endl;
  BasicFunction my_binomial(binomial);
  cout << "The functor at zero is: " << my_binomial(0) << "." << endl;
  cout << "Finding the root." << endl;
  double root = bisection_root_finder(my_binomial, XMIN,XMAX);
  cout << "The root is: " << root << endl;
  cout << "The functor at the root is: " << my_binomial(root) << endl;
  cout << "The function was called "
	    << my_binomial.get_times_called() << " times." << endl;

  return 0;
}
