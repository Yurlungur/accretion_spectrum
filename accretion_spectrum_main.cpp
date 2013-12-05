// accretion_spectum_main.cpp

// Time-stamp: <2013-12-05 17:00:42 (jonah)>
// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

// This is the main loop for the accretion spectrum solver.
// ----------------------------------------------------------------------

// Includes
#include "accretion_spectrum.hpp"
#include <fstream>
using std::cout;
using std::endl;
using std::abs;
// ----------------------------------------------------------------------

int main() {
  // Helpful message
  cout << "Solving for the temperature of an accretion disk around a\n"
       << "pseudo-Newtonian stellar-mass black hole." << endl;
  cout << "Parameters are:\n"
       << "\tAdiabatic index: " << GAMMA << "\n"
       << "\talpha: " << ALPHA << "\n"
       << "\tf: " << F << "\n"
       << "\tspin a: " << A << "\n"
       << "\tNumber of solar masses: " << NUM_SOLAR_MASSES << "\n"
       << "\tFraction eddington accretion: " << NUM_EDDINGTON_RATES << "\n"
       << "\tOutput file name: " << OUTPUT_FILE_NAME
       << endl;

  // Initial guesses for max and min values of Rs
  double r_s_min = 2.5*R_G;
  double r_s_max = 2.5*1E2*R_G;

  // Helpful message
  cout << "Bounds on sonic point:\n"
       << "\t[" << r_s_min << ", " << r_s_max << "]"
       << endl;

  cout << "Working..." << endl;

  // Important data for the while loop.
  bool angular_momentum_condition_satisfied = false;
  double Rs;
  bool all_okay;

  // Make an integrator
  Integrator rk_integrator;
  RKF45 outward_integrator;

  // Now we want to find roots until the angular momentum condition is
  // satisfied
  while ( !angular_momentum_condition_satisfied ) {
    Rs = bisection_root_finder(rk_integrator, r_s_min, r_s_max);
    angular_momentum_condition_satisfied
      = abs(angular_derivative_condition(rk_integrator.integrator))
      < ROOT_ERROR;
    if ( !angular_momentum_condition_satisfied ) {
      r_s_max = Rs;
    }
  }
  // Helpful message.
  cout << "Found sonic point: " << Rs << endl;

  // Now integrate all the way inward.
  cout << "Integrating inward to event horizon..." << endl;
  finish_integration(rk_integrator);

  // Once we know both conditions are satisfied, we want to save the
  // data to a file.
  cout << "Saving data..." << endl;
  print_inward_integrator_to_file(rk_integrator, OUTPUT_FILE_NAME);

  // And we want to check to make sure this is really the integrator
  // we want.
  cout << "As a last check, ensuring the data is consistent." << endl;
  all_okay = outward_integration_okay(rk_integrator,outward_integrator);
  if (all_okay) {
    cout << "Great! Everything worked! Have a nice day." << endl;
  } else {
    cout << "Oh dear. Something went wrong. Dumping integration data."
	 << endl;
    std::ofstream file;
    file.open(INTEGRATION_DUMP);
    file << outward_integrator << endl;
    file.close();
  }
  return 0;
}
