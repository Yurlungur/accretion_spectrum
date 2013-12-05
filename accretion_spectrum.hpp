// accretion_spectrum.hpp

// Time-stamp: <2013-12-05 17:02:53 (jonah)>
// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

// This is the prototype of the simulation of an axisymmetric
// accretion disk.

// Include guard
#pragma once

// Includes
#include "bisection.hpp"

// Constants
// ----------------------------------------------------------------------
// Parameters we might want to change
const double GAMMA = 1.5; // The adiabatic index of our gas
const double ALPHA = 0.1; // The Shakura-Sunyaev alpha
const double F = 3.5E-5; // (1-F) is the radiative efficiency of the disk.
const double A = 0; // The unit-less spin parameter for black hole.
const double NUM_SOLAR_MASSES = 10; // The number of masses our black hole is
const double NUM_EDDINGTON_RATES = 0.1; // How much of the eddington
					// accretion rate we want our
					// accretion rate to be.
const double ROUT_OVER_RS = 1E5; // how big rout is compared to Rs.
const int NUM_VARIABLES = 3; // The number of variables we
			     // have: v_R,omega,c2s
const char OUTPUT_FILE_NAME [] = "schwarzschild_bh.dat";
const char INTEGRATION_DUMP [] = "integration_dump.dat";

// Fundamental constants. Uses mks
const double C = 299792458; // m/s the speed of light
const double G = 6.67384E-11; // m^3 /kg s^2. Newton's constant
const double M_SUN = 1.9891E30; // kg. 1 solar mass
const double M_PROTON = 1.672621777E-27; // kg. The mass of a proton
const double SIGMA_T = 6.6524587E-29; // m^2. The Thomson scattering
				      // cross section
const double PLANCK_CONSTANT = 6.62606957E-34; // Joule-seconds. Planck's
					       // constant
const double SIGMA_B = 5.670400E-8; // J s^{-1} m^{-2} K^{-4}. The
				    // Stephan-Boltzmann Constant.
const double BOLTZMANN_CONSTANT = 1.38046488E-23; // J/K. Boltzmann constant

// Constants related to our simulation
const double M = NUM_SOLAR_MASSES*M_SUN; // the mass of our black hole
// The eddington accretion limit
const double M_DOT_EDDINGTON = (4*M_PI*G*M*M_PROTON)/((1-F)*C*SIGMA_T);
// Our accretion rate
const double M_DOT = -1 * NUM_EDDINGTON_RATES * M_DOT_EDDINGTON;
const double R_G = (G*M)/(C*C); // Twice the Schwarzschild radius of a black hole

// Self-similarity constants
const double EPSILON_PRIME = ((5.0/3.0) - GAMMA)/(F * (GAMMA - 1));
const double C2_0 = 2/(5 + 2*EPSILON_PRIME + (ALPHA * ALPHA/EPSILON_PRIME));
const double V_0 = -ALPHA*sqrt(C2_0/EPSILON_PRIME);
const double OMEGA_0 = sqrt(C2_0 * EPSILON_PRIME);

// Constants relating to the Runge-Kutta integrator
const double MIN_NUMBER_POINTS = 1000; // Max number of points sampled
				       // with the Runge-Kutta
				       // integrator
const double INTEGRATOR_RELATIVE_ERROR_FACTOR = 1E-5; // what fraction
						      // of the
						      // solution is
						      // allowed to
						      // have what
						      // error
const double SCHWARZSCHILD_CLOSENESS = 0.9; // How close to the
					    // Schwarzschild radius we
					    // integrate.

// ----------------------------------------------------------------------


// For convenience, we overload the unary minus-sign operator for
// dVectors.
// ----------------------------------------------------------------------
dVector operator-(const dVector& in); 
// ----------------------------------------------------------------------


// Helper functions
// ----------------------------------------------------------------------
// The Keplerian angular velocity of the black hole
double omega_k(double R);
// The derivative of the Keplerian angular velocity of the black hole
double omega_k_prime(double R);
// R_out as a function of Rs
double get_R_out(double Rs);
// The new time variable, tau. Depends on R and Rs.
double get_tau(double R,double Rs);
// R as a function of tau and Rs.
double get_R(double tau, double Rs);
// Get H from the other variables
double get_H(double R, double c2s);
// Get rho from the other variables
double get_rho(double R, double v_R, double c2s);
// ----------------------------------------------------------------------


// The self-similar boundary data
// ----------------------------------------------------------------------
double v_R_ss(double R_out); // radial velocity
double c2s_ss(double R_out); // square speed of sound 
double omega_ss(double R_out); // angular velocity
// ----------------------------------------------------------------------


// The sonic-point boundary data
// ----------------------------------------------------------------------
double speed_condition(double v_R, double c2s);
double angular_derivative_condition(double R, double v_R, double c2s,
				    double omega, double omega_prime);
double speed_condition(const RKF45& integrator);
double angular_derivative_condition(const RKF45& integrator);
// ----------------------------------------------------------------------


// The equations for the derivatives
// ----------------------------------------------------------------------
// This method is just a helper function this denominator is the same
// for every one of the odes
double get_ode_denominator(double R, double v_R, double c2s, double omega);
double get_v_R_prime(double R, double v_R, double c2s, double omega);
double get_c2s_prime(double R, double v_R, double c2s, double omega);
double get_omega_prime(double R, double v_R, double c2s, double omega);
// ----------------------------------------------------------------------


// The vector y = [v_R,c2s,omega] gives us the function to iterate,
// y' = f(R,y).
// We also want another function, though, z' = -f(tau,z). This is to
// integrate from large R to small. Both methods are defined here.
// ----------------------------------------------------------------------
dVector get_y_prime(double R, const dVector& y);
// Rs is the sonic point. Needed to calculate tau. It's stored in a
// zero-dimensional dVector for convenience. This is how the
// integrator expects it.
dVector get_z_prime(double tau, const dVector& z, const dVector& Rs_vector);
// ----------------------------------------------------------------------


// In the end, we want to calculate a spectrum, so the following are
// useful calculations to do in C++.
// ----------------------------------------------------------------------
double get_outward_torque(double R, double v_R, double c2s,
			  double omega_prime);
double get_dissipated_energy(double R, double v_R, double c2s,
			     double omega_prime);
double get_temperature(double R, double v_R, double c2s, double omega_prime);
// nu is the frequency of the light. In Hz. 
double get_panck_intensity(double nu, double R, double v_R, double c2s,
			   double omega_prime);
// ----------------------------------------------------------------------


// Integrator subclass. Specifically for the shooting method.
// ---------------------------------------------------------------------- 
class Integrator: public ScalarFunction {
public:
  // Constructor
  Integrator();
  // Destructor
  ~Integrator() {}

  // The input/output. x is the guess for the sonic point, Rs. The
  // integrator handles the rest.
  double operator()(double x);

  // The integrator used. Uninitialized.
  RKF45 integrator;

protected:
  double max_t;

  // Initialize an integrater with the appropriate initial data for a
  // guess of Rs.
  void initialize_integrator(double Rs);
  // Evaluate the speed condition on the integrator.
  double evaluate_speed_condition();

};
// ----------------------------------------------------------------------


// Other integrator stuff
// ----------------------------------------------------------------------
// Once we solve the BVP integrating inwards, we want to double check
// by integrating outwards again this method double checks. It also
// spits out an RKF45 integrator that has completed the integration.
bool outward_integration_okay(const Integrator& inward_Integrator,
			      RKF45& outward_integrator);

// Finsh integrating inwards to the Schwarzschild radius.
void finish_integration(Integrator& inward_Integrator);
// ----------------------------------------------------------------------


// Output
// ----------------------------------------------------------------------
// Prints the integrator to a stream with all the useful
// information. Assumes the integrator is inward and uses the variable
// tau not R.
void print_inward_integrator(const Integrator& integrator, std::ostream& s);
// Prints the same integrator to a file
void print_inward_integrator_to_file(const Integrator& integrator,
				     const char* filename);
// ----------------------------------------------------------------------
