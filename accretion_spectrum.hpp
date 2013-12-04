// accretion_spectrum.hpp

// Time-stamp: <2013-12-04 18:55:24 (jonah)>
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
const char OUTPUT_FILE_NAME [] = "schwarzschild_bh.dat";

// Fundamental constants. Uses mks
const double C = 299792458; // m/s the speed of light
const double G = 6.67384E-11; // m^3 /kg s^2. Newton's constant
const double M_SUN = 1.9891E30; // kg. 1 solar mass
const double M_PROTON = 1.672621777E-27; // kg. The mass of a proton
const double SIGMA_T = 6.6524587E-29; // m^2. The Thomson scattering
				      // cross section
const double PLANCK_CONSTANT = 6.62606957E-34; // Joule-seconds. Planck's constant

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
// ----------------------------------------------------------------------
