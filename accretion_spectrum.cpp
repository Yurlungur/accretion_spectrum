// accretion_spectrum.cpp

// Time-stamp: <2013-12-04 18:53:29 (jonah)>
// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

// This is the implementation of an axisymmetric accretion disk.
// ----------------------------------------------------------------------

// Includes
#include "accretion_spectrum.hpp"


// Helper functions
// ----------------------------------------------------------------------
double omega_k(double R) {
  double r = R/R_G;
  return (pow(C,3)/(G*M))*((r*r - 2*A*sqrt(r) + A*A)/(r*r*(r-2) + A));
}
double omega_k_prime(double R) {
  double r = R/R_G;
  double denominator = 2*sqrt(r)*M*pow((R-2*R_G)*sqrt(r) + A*R_G,2)*G*pow(R,3);
  double numerator = pow(C,3)*(-3*pow(R,4)
			       + 12*pow(r,1.5)*pow(R_G,3)*A*R
			       - 16*pow(r,1.5)*pow(R_G,4)*A
			       + 16*pow(R_G,3)*R*A*A
			       - 7*A*A*R_G*R_G*R*R
			       - 4*sqrt(r)*pow(R_G,4)*pow(A,3)
			       + 2*pow(R,3)*R_G);
  return numerator/denominator;
}

double get_R_out(double Rs) {
  return ROUT_OVER_RS * Rs;
}
double get_tau(double R, double Rs) {
  return (R - get_R_out(Rs))/(Rs - get_R_out(Rs));
}
double get_R(double tau, double Rs) {
  return get_R_out(Rs) + (Rs - get_R_out(Rs))*tau;
}
double get_H(double R, double c2s) {
  return c2s/omega_k(R);
}
double get_rho(double R, double v_R, double c2s) {
  return -M_DOT/(4*M_PI*v_R*R*get_H(R,c2s));
}
// ----------------------------------------------------------------------


// The self-similar solutions
// ----------------------------------------------------------------------
double v_R_ss(double R_out) {
  return V_0*sqrt(G*M/R_out);
}
double c2s_ss(double R_out) {
  return C2_0*G*M/R_out;
}
double omega_ss(double R_out) {
  return OMEGA_0 * omega_k(R_out);
}
// ----------------------------------------------------------------------


// The sonic-point boundary data
// ----------------------------------------------------------------------
double speed_condition(double v_R, double c2s) {
  return v_R*v_R - (2*GAMMA*c2s)/(GAMMA+1);
}
double angular_derivative_condition(double R, double v_R, double c2s,
				    double omega, double omega_prime) {
  double omega_kp = omega_k_prime(R);
  double term1,term2,term3;
  term1 = (pow(omega_k(R),2) - pow(omega,2))*R;
  term2 = ((c2s*2*GAMMA)/(GAMMA+1))*((1.0/R)-(1.0/omega_k(R))*omega_kp);
  term3 = c2s*((GAMMA-1)/(GAMMA+1))*((F*ALPHA*R)/v_R)*omega_prime;
  return term1 - term2 - term3;
}
// ----------------------------------------------------------------------
