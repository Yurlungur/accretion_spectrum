// accretion_spectrum_dry_run.cpp

// Time-stamp: <2013-12-08 15:27:39 (jonah)>
// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

// This is a sort of "test run" for the accretion spectrum solver. We
// integrate inward for a "reasonable guess" for the sonic
// point. Additionally, we evaluate both constraints for a range of
// guesses of the sonic radius.

// ----------------------------------------------------------------------

// Includes
#include "accretion_spectrum.hpp"
#include <fstream>
using std::cout;
using std::endl;
using std::abs;
// ----------------------------------------------------------------------


// Constants
// ----------------------------------------------------------------------
const int NUM_POINTS = 1000; // number of points in constraint plot
// Name of integration test file
const char INTEGRATION_TEST [] = "integration_test.dat";
// Name of constraint plot file
const char CONSTRAINT_PLOT [] = "constraints_plot.dat";
// ----------------------------------------------------------------------


// main function
// ----------------------------------------------------------------------
int main() {
  double rs_guess = 2.3*R_G;
  double r_out_guess = get_R_out(rs_guess);
  double rs_min = 2.4*R_G;
  double rs_max = NUM_POINTS*rs_min;
  double rs_evals [NUM_POINTS];
  double speed_test_results [NUM_POINTS];
  double angular_test_results [NUM_POINTS];
  double my_v_R_ss,my_c2s_ss,my_omega_ss;
  double my_v_R_prime,my_c2s_prime,my_omega_prime;
  dVector initial_data;
  dVector optional_args;
  Integrator rk_integrator;
  RKF45 test_integrator;
  std::ofstream file;

  // Helpful message
  cout << "Solving for the temperature of an accretion disk around a\n"
       << "pseudo-Newtonian stellar-mass black hole." << endl;
  cout << "Parameters are:\n"
       << "\tAdiabatic index: " << GAMMA << "\n"
       << "\talpha: " << ALPHA << "\n"
       << "\tf: " << F << "\n"
       << "\tspin a: " << A << "\n"
       << "\tNumber of solar masses: " << NUM_SOLAR_MASSES << "\n"
       << "\tFraction eddington accretion: " << NUM_EDDINGTON_RATES
       << endl;

  cout << "Integrating inward and guessing that the sonic radius is "
       << rs_guess << "." << endl;
  rk_integrator(rs_guess);
  cout << "Printing results to " << INTEGRATION_TEST << "." << endl;
  print_inward_integrator_to_file(rk_integrator,INTEGRATION_TEST);

  cout << "Trying to integrate inward with just the RK45 class."
       << endl;
  test_integrator.set_f(get_z_prime);
  test_integrator.set_relative_error_factor(INTEGRATOR_RELATIVE_ERROR_FACTOR);
  initial_data.push_back(v_R_ss(r_out_guess));
  initial_data.push_back(c2s_ss(r_out_guess));
  initial_data.push_back(omega_ss(r_out_guess));
  test_integrator.set_y0(initial_data);
  optional_args.push_back(rs_guess);
  test_integrator.set_t0(0);
  test_integrator.set_max_dt(1.0/MIN_NUMBER_POINTS);
  test_integrator.set_optional_args(optional_args);
  test_integrator.integrate(1);
  cout << "Feeding to file " << INTEGRATION_DUMP << "." << endl;
  file.open(INTEGRATION_DUMP);
  file << test_integrator;
  file.close();

  cout << "Finding the values of the constraints at the sonic radius\n"
       << "for various guesses of the sonic radius\n"
       << "between " << rs_min << " and " << rs_max << "."
       << endl;
  for (int i = 0; i < NUM_POINTS; i++) {
    rs_evals[i] = (i+1)*rs_min;
    speed_test_results[i] = rk_integrator(rs_evals[i]);
    angular_test_results[i] = rk_integrator.evaluate_angular_condition();

    my_v_R_ss = v_R_ss(rs_evals[i]);
    my_c2s_ss = c2s_ss(rs_evals[i]);
    my_omega_ss = omega_ss(rs_evals[i]);
    cout << "Self similar data:\n"
	 << "\tv_R = " << my_v_R_ss << "\n"
	 << "\tc2s = " << my_c2s_ss << "\n"
	 << "\tomega = " << my_omega_ss << "." << endl;
    my_v_R_prime = get_v_R_prime(rs_evals[i],my_v_R_ss,
				 my_c2s_ss,my_omega_ss);
    my_c2s_prime = get_c2s_prime(rs_evals[i],my_v_R_ss,
				 my_c2s_ss,my_omega_ss);
    my_omega_prime = get_omega_prime(rs_evals[i],my_v_R_ss,
				     my_c2s_ss,my_omega_ss);
    cout << "Derivatives assuming self-similar data values:\n"
	 << "\tv_R = " << my_v_R_prime << "\n"
	 << "\tc2s = " << my_c2s_prime << "\n"
	 << "\tomega = " << my_omega_prime << "." << endl;
  }

  cout << "Printing to file " << CONSTRAINT_PLOT << "." << endl;
  file.open(CONSTRAINT_PLOT);
  file << "# Guess for Rs\tRadial constraint\tAngular constraint" << endl;
  for (int i = 0; i < NUM_POINTS; i++) {
    file << rs_evals[i] << "\t"
	 << speed_test_results[i] << "\t"
	 << angular_test_results[i] << endl;
  }
  file.close();

  return 0;
}
// ----------------------------------------------------------------------
