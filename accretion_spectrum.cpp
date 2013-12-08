// accretion_spectrum.cpp

// Time-stamp: <2013-12-08 01:31:24 (jonah)>
// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

// This is the implementation of an axisymmetric accretion disk.
// ----------------------------------------------------------------------

// Includes
#include "accretion_spectrum.hpp"
#include <cassert>
#include <fstream>
using std::abs;
// ----------------------------------------------------------------------


// For convenience, we overload the unary minus-sign operator for
// dVectors.
// ----------------------------------------------------------------------
dVector operator-(const dVector& in) {
  dVector out(in.size());
  for (int i = 0; i < (int)in.size(); i++) {
    out[i] = -in[i];
  }
  return out;
}
// ----------------------------------------------------------------------



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
  return V_0*sqrt((G*M)/R_out);
}
double c2s_ss(double R_out) {
  return (C2_0*(G*M))/R_out;
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
double speed_condition(const RKF45& integrator) {
  dVector y_vec = integrator.get_y();
  double v_R = y_vec[0];
  double c2s = y_vec[1];
  return speed_condition(v_R,c2s);
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
double angular_derivative_condition(const RKF45& integrator) {
  dVector y_vec = integrator.get_y();
  double v_R = y_vec[0];
  double c2s = y_vec[2];
  double omega = y_vec[3];
  double Rs = integrator.get_optional_args()[0];
  double tau = integrator.get_t();
  double R = get_R(tau,Rs);
  double omega_prime = get_omega_prime(R,v_R,c2s,omega);
  return angular_derivative_condition(R,v_R,c2s,omega,omega_prime);
}
// ----------------------------------------------------------------------



// The equations for the derivatives
// ----------------------------------------------------------------------
void set_psi_chi(double& psi, double& chi, double c2s) {
  psi = ALPHA*c2s;
  chi = F*psi*psi;
}

double get_L(double R, double omega) {
    return R*R*(omega*omega - pow(omega_k(R),2));
}

double get_gamma(double R, double v_R, double c2s, double omega) {
  return 0.5*R*omega_k(R)*(pow(v_R,4)*(GAMMA + 1)
			   + 2*v_R*v_R*c2s*(F*pow(ALPHA,2)*(GAMMA - 1)
					    + GAMMA)
			   - F*(GAMMA-1)*pow(ALPHA*c2s,2));
}

double get_v_R_prime(double R, double v_R, double c2s, double omega) {
  double numerator, denominator,psi,chi,sigma,L;
  denominator = 2*get_gamma(R,v_R,c2s,omega);
  set_psi_chi(psi,chi,c2s);
  L = get_L(R,omega);
  sigma = -v_R*v_R*GAMMA + F*ALPHA*R*v_R*omega*(GAMMA-1)
    + F*ALPHA*ALPHA*L*(GAMMA-1);
  numerator = v_R*(omega_k(R)*(3*(GAMMA-1)*chi + 2*c2s*sigma
			       - v_R*v_R*L*(GAMMA+1))
		   -2*R*c2s*omega_k_prime(R)*(F*ALPHA*ALPHA*c2s*(GAMMA-1)
					      - v_R*v_R*GAMMA));
  return numerator/denominator;
}

double get_c2s_prime(double R, double v_R, double c2s, double omega) {
  double numerator,denominator,psi,chi,L,xi,prefactor;
  denominator = get_gamma(R,v_R,c2s,omega);
  set_psi_chi(psi,chi,c2s);
  L = get_L(R,omega);
  xi = pow(v_R,2) - 2*F*ALPHA*R*v_R*omega + L;
  prefactor = (GAMMA-1)*c2s;
  numerator = prefactor*(omega_k(R)*(2*chi +
				     F*psi*(ALPHA*L
					    + 2*R*v_R*omega
					    - ALPHA*v_R*v_R)
				     + v_R*v_R*xi)
			 - R*omega_k_prime(R)*(chi - pow(v_R,4)));
  return numerator/denominator;
}

double get_omega_prime(double R, double v_R, double c2s, double omega) {
  double numerator,denominator,psi,chi,L,lambda;
  denominator = 2*R*get_gamma(R,v_R,c2s,omega);
  set_psi_chi(psi,chi,c2s);
  L = get_L(R,omega);
  lambda = ALPHA*v_R*v_R*(GAMMA - 3)
    + 4*R*omega*v_R*GAMMA
    + ALPHA*L*(3*GAMMA - 1);
  numerator = -v_R*(omega_k(R)*(4*ALPHA*c2s*c2s*GAMMA
				+ c2s*lambda
				- 2*R*pow(v_R,3)*omega*(GAMMA+1))
		   -2*ALPHA*R*c2s*omega_k_prime(R)*(c2s*GAMMA
						    +v_R*v_R*(GAMMA-1)));
  return numerator/denominator;
}
// ----------------------------------------------------------------------


// The vector y = [v_R,c2s,omega] gives us the function to iterate,
// y' = f(R,y).
// We also want another function, though, z' = -f(tau,z). This is to
// integrate from large R to small. Both methods are defined here.
// ----------------------------------------------------------------------
dVector get_y_prime(double R, const dVector& y) {
  assert ( (int)y.size() == NUM_VARIABLES
	   && "There are three variables we're integrating over." );
  double v_R,c2s,omega;
  dVector y_prime(NUM_VARIABLES);
  // convenience bindings
  v_R = y[0];
  c2s = y[1];
  omega = y[2];
  // Get the primes
  y_prime[0] = get_v_R_prime(R,v_R,c2s,omega);
  y_prime[1] = get_c2s_prime(R,v_R,c2s,omega);
  y_prime[2] = get_omega_prime(R,v_R,c2s,omega);
  
  return y_prime;
}

dVector get_z_prime(double tau, const dVector& z, const dVector& Rs_vector) {
  assert ( Rs_vector.size() == 1 && "There's only one Rs!" );
  double Rs = Rs_vector[0];
  double R = get_R(tau,Rs);
  return -get_y_prime(R,z);
}
// ----------------------------------------------------------------------


// In the end, we want to calculate a spectrum, so the following are
// useful calculations to do in C++.
// ----------------------------------------------------------------------
double get_outward_torque(double R, double v_R, double c2s,
			  double omega_prime) {
  double rho = get_rho(R,v_R,c2s);
  double H = get_H(R, c2s);
  return - ALPHA*rho*sqrt(c2s)*pow(H*R,2)*omega_prime;
}
double get_dissipated_energy(double R, double v_R, double c2s,
			     double omega_prime) {
  return get_outward_torque(R,v_R,c2s,omega_prime)*omega_prime/(4*M_PI*R);
}
double get_temperature(double R, double v_R, double c2s, double omega_prime) {
  double j = -(1-F)*get_dissipated_energy(R,v_R,c2s,omega_prime);
  /*
  if (DEBUGGING) {
    std::cout << "\tj = " << j << std::endl;;
    std::cout << "\tT = " << pow(j/SIGMA_B,0.25) << std::endl;
  }
  */
  return pow(j/SIGMA_B,0.25);
}
// nu is the frequency of the light. In Hz. 
double get_panck_intensity(double nu, double R, double v_R, double c2s,
			   double omega_prime) {
  double prefactor = 2*PLANCK_CONSTANT*pow(nu,3)*C;
  double T = get_temperature(R,v_R,c2s,omega_prime);
  double exponent = (PLANCK_CONSTANT*nu)/(BOLTZMANN_CONSTANT*T);
  return prefactor/(exp(exponent) - 1);
}
// ----------------------------------------------------------------------


// Integrator subclass. Specifically for the shooting method.
// ----------------------------------------------------------------------
Integrator::Integrator() {
  integrator = RKF45(get_z_prime);
  max_t = 1;
  integrator.set_relative_error_factor(INTEGRATOR_RELATIVE_ERROR_FACTOR);
}
// The input/output. x is the guess for the sonic point, Rs. The
// integrator handles the rest.
double Integrator::operator()(double x) {
  initialize_integrator(x);
  integrator.integrate(max_t);
  double output = evaluate_speed_condition();
  dVector data;

  if (DEBUGGING) {
    data = integrator.get_y();
    std::cout << "Data: "
	      << "\tv_R: " << data[0] << "\n"
	      << "\tc_s: " << data[1] << "\n"
	      << "\tomega: " << data[2] << "." << std::endl;
  }

  return output;
}

void Integrator::initialize_integrator(double Rs) {
  // Initial data
  double R_out = get_R_out(Rs);
  double v_R_initial = v_R_ss(R_out);
  double c2s_initial = c2s_ss(R_out);
  double omega_initial = omega_ss(R_out);
  dVector initial_data(NUM_VARIABLES);
  initial_data[0] = v_R_initial;
  initial_data[1] = c2s_initial;
  initial_data[2] = omega_initial;
  integrator.set_y0(initial_data);

  // Integrates from t = tau = 0 to t = tau = 1.
  // This goes from R_out to Rs.x
  dVector optional_args(1);
  optional_args[0] = Rs;
  integrator.set_t0(0);
  integrator.set_max_dt(1.0/MIN_NUMBER_POINTS);
  integrator.set_optional_args(optional_args);
  // Resets the integrator, deleting all known information about time
  // evolution
  integrator.reset();
}

double Integrator::evaluate_speed_condition() {
  return speed_condition(integrator);
}

double Integrator::evaluate_angular_condition() {
  return angular_derivative_condition(integrator);
}
//----------------------------------------------------------------------


// Other integrator stuff
// ----------------------------------------------------------------------
bool outward_integration_okay(const Integrator& inward_Integrator,
			      RKF45& outward_integrator) {
  // Let's get some variables out of the way
  double Rs = inward_Integrator.integrator.get_optional_args()[0];
  double R_out = get_R_out(Rs);
  double v_R_out_condition = v_R_ss(R_out);
  double c2s_out_condition = c2s_ss(R_out);
  double omega_out_condition = omega_ss(R_out);
  // Set the integrator
  outward_integrator.set_f(get_y_prime);
  outward_integrator.set_y0(inward_Integrator.integrator.get_y());
  outward_integrator.set_t0(Rs);
  outward_integrator.set_relative_error_factor(INTEGRATOR_RELATIVE_ERROR_FACTOR);
  outward_integrator.set_max_dt((R_out-Rs)/MIN_NUMBER_POINTS);
  outward_integrator.reset();
  // Integrate, baby!
  outward_integrator.integrate(R_out);
  // Test whether or not we're okay
  dVector final_y = outward_integrator.get_y();
  double v_R_out = final_y[0];
  double c2s_out = final_y[1];
  double omega_out = final_y[2];
  return ( abs(v_R_out - v_R_out_condition) < ROOT_ERROR
	   && abs(c2s_out - c2s_out_condition) < ROOT_ERROR
	   && abs(omega_out - omega_out_condition) < ROOT_ERROR );
}

// Finsh integrating inwards to the Schwarzschild radius.
void finish_integration(Integrator& inward_Integrator) {
  double Rs = inward_Integrator.integrator.get_optional_args()[0];
  double final_tau = SCHWARZSCHILD_CLOSENESS*get_tau(R_G,Rs);
  inward_Integrator.integrator.integrate(final_tau);
}
// ----------------------------------------------------------------------


// Output
// ----------------------------------------------------------------------
void print_inward_integrator(const Integrator& integrator,
			     std::ostream& out) {
  double Rs = integrator.integrator.get_optional_args()[0];
  dVector y;
  double R,v_R,c2s,omega,rho,omega_prime,H,T;
  out << "# R v_R c2s omega rho H omega_prime T" << std::endl;
  for (int n = 0; n < integrator.integrator.steps(); n++) {
    y = integrator.integrator.get_y(n);
    R = get_R(integrator.integrator.get_t(n),Rs);
    v_R = y[0];
    c2s = y[1];
    omega = y[2];
    omega_prime = get_omega_prime(R,v_R,c2s,omega);
    rho = get_rho(R,v_R,c2s);
    H = get_H(R,c2s);
    T = get_temperature(R,v_R,c2s,omega_prime);
    out << R << " " << v_R << " " << c2s << " " << omega << " "
	<< rho << " " << H << " " << omega_prime << " " << T
	<< std::endl;
  }
}

void print_inward_integrator_to_file(const Integrator& integrator,
				     const char* filename) {
  std::ofstream file;
  file.open(filename);
  print_inward_integrator(integrator,file);
  file.close();
}
// ----------------------------------------------------------------------
