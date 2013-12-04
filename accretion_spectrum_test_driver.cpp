// accretion_spectrum_test_driver.cpp

// Time-stamp: <2013-12-04 18:17:53 (jonah)>
// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

// This is the unit testing toolbox for the accretion-spectrum
// package.

#include "accretion_spectrum.hpp"
using namespace std;

const double R_TEST = 10*R_G;

int main() {
  cout << "Testing accretion-spectrum functions\n"
       << "---------------------------------------------------------\n"
       << endl;

  cout << "First, test the Keplerian angular velocity methods.\n"
       << "---------------------------------------------------------"
       << endl;

  cout << "Using R_TEST = " << R_TEST << ",\n"
       << "the Keplerian angular velocity and its derivative are:\n"
       << "\t Omega_k = " << omega_k(R_TEST) << "\n"
       << "\t Omega_k_prime = " << omega_k_prime(R_TEST)
       << endl;

  return 0;
}
