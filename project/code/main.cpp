#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "constants.h"
#include "../../mclib/mclib.h"
#include "radiation.h"
#include "utilities.h"
#include "configuration.h" // Include last to allow redefinitions.

void rhs(double t, double const * Y_0, double * R);

int main(int argc, char * argv[]) {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	
	// Configure vertical coordinates.
	double z_TOA, dz; // / m
	double * z; // / m
	double * P; // / Pa
	double * delta, * sigma;
	z_TOA = get_altitude(global_P_TOA);
	dz = (z_TOA - global_z_g) / global_N;
	z = new double[global_N + 1];
	delta = new double[global_N + 1];
	P = new double[global_N + 1];
	sigma = new double[global_N + 1];
	for (int i = 0; i <= global_N; i++) {
		z[i] = z_TOA - i * dz;
		delta[i] = get_optical_depth_z(z[i], global_P_TOA);
		P[i] = get_pressure(z[i]);
		sigma[i] = get_sigma(P[i], global_P_TOA);
	}
	
	
	
	/* Analytical solution in radiative equilibrium */
	
	// Prepare variables.
	double * T, * theta, * E_U, * E_D;
	double T_0; // / K
	T = new double[global_N + 1];
	theta = new double[global_N + 1];
	E_U = new double[global_N + 1];
	E_D = new double[global_N + 1];
	T_0 = pow(global_S_t / global_sigma, 0.25);
	
	// Prepare output files.
	ofstream file_temperature, file_irradiance;
	char fn_temperature_analytical[] = DIR_DATA "/temperature_analytical.dat";
	char fn_irradiance_analytical[] = DIR_DATA "/irradiance_analytical.dat";
	file_temperature << fixed << setprecision(N_PRECISION);
	file_temperature.open(fn_temperature_analytical);
	file_temperature << "#z/z_0 T delta P sigma theta" << endl;
	file_temperature << "#'1' 'K' '1' 'Pa' '1' 'K'" << endl;
	file_irradiance << fixed << setprecision(N_PRECISION);
	file_irradiance.open(fn_irradiance_analytical);
	file_irradiance << "#z/z_0 E_U/S_t E_D/S_t delta P sigma" << endl;
	file_irradiance << "#'1' '1' '1' '1' 'Pa' '1'" << endl;
	
	// Plot analytical solutions.
	for (int i = 0; i <= global_N; i++) {
		T[i] = temperature_norm(delta[i]);
		theta[i] = get_theta(T[i], P[i]);
		file_temperature << z[i] / global_z_0 << ' ' << T_0 * T[i] << ' ' << delta[i] << ' ' << P[i] << ' ' << sigma[i] << ' ' << T_0 * theta[i] << '\n';
		E_U[i] = irradiance_upward_norm(delta[i]);
		E_D[i] = irradiance_downward_norm(delta[i]);
		file_irradiance << z[i] / global_z_0 << ' ' << E_U[i] << ' ' << E_D[i] << ' ' << delta[i] << ' ' << P[i] << ' ' << sigma[i] << '\n';
	}
	file_temperature.close();
	cout << "Analytical solution temperature profile in radiative equilibrium calculated, values are stored in file " << fn_temperature_analytical << endl;
	file_irradiance.close();
	cout << "Analytical solution irradiances in radiative equilibrium calculated, values are stored in file " << fn_irradiance_analytical << endl;
	
	
	
	/* Numerical solution in radiative equilibrium */
	
	// Prepare variables and set initial conditions.
	double Y[3];
	double T_tmp, theta_tmp;
	Y[0] = 0.5;
	Y[1] = 1.0;
	Y[2] = 0.0;
	T_tmp = pow(Y[0], 0.25);
	theta_tmp = get_theta(T_tmp, P[0]);
	
	// Prepare output files.
	char fn_temperature_numerical[] = DIR_DATA "/temperature_numerical.dat";
	char fn_irradiance_numerical[] = DIR_DATA "/irradiance_numerical.dat";
	ofstream file_errors;
	char fn_errors[] = DIR_DATA "/errors.dat";
	file_temperature << fixed << setprecision(N_PRECISION);
	file_temperature.open(fn_temperature_numerical);
	file_temperature << "#z/z_0 T delta P sigma theta" << endl;
	file_temperature << "#'1' 'K' '1' 'Pa' '1' 'K'" << endl;
	file_irradiance << fixed << setprecision(N_PRECISION);
	file_irradiance.open(fn_irradiance_numerical);
	file_irradiance << "#z/z_0 E_U/S_t E_D/S_t delta P sigma" << endl;
	file_irradiance << "#'1' '1' '1' '1' 'Pa' '1'" << endl;
	file_errors << scientific << setprecision(N_PRECISION);
	file_errors.open(fn_errors);
	file_errors << "#z/z_0 T/T_0 theta/T_0 E_U/S_t E_D/S_t delta P sigma" << endl;
	file_errors << "#'1' '1' '1' '1' '1' '1' 'Pa' '1'" << endl;
	
	// Integrate equations.
	file_temperature << z[0] / global_z_0 << ' ' << T_0 * T_tmp << ' ' << delta[0] << ' ' << P[0] << ' ' << sigma[0] << ' ' << T_0 * theta[0] << '\n';
	file_irradiance << z[0] / global_z_0 << ' ' << Y[1] << ' ' << Y[2] << ' ' << delta[0] << ' ' << P[0] << ' ' << sigma[0] << '\n';
	file_errors << z[0] / global_z_0 << ' ' << fabs(T_tmp - T[0]) << ' ' << fabs(theta_tmp - theta[0]) << ' ' << fabs(Y[1] - E_U[0]) << ' ' << fabs(Y[2] - E_D[0]) << ' ' << delta[0] << ' ' << P[0] << ' ' << sigma[0] << '\n';
	for (int i = 1; i <= global_N; i++) {
		rungekutta4(delta[i], delta[i] - delta[i-1], Y, rhs, 3);
		T_tmp = pow(Y[0], 0.25);
		theta_tmp = get_theta(T_tmp, P[i]);
		file_temperature << z[i] / global_z_0 << ' ' << T_0 * T_tmp << ' ' << delta[i] << ' ' << P[i] << ' ' << sigma[i] << ' ' << T_0 * theta[i] << '\n';
		file_irradiance << z[i] / global_z_0 << ' ' << Y[1] << ' ' << Y[2] << ' ' << delta[i] << ' ' << P[i] << ' ' << sigma[i] << '\n';
		file_errors << z[i] / global_z_0 << ' ' << fabs(T_tmp - T[i]) << ' ' << fabs(theta_tmp - theta[i]) << ' ' << fabs(Y[1] - E_U[i]) << ' ' << fabs(Y[2] - E_D[i]) << ' ' << delta[i] << ' ' << P[i] << ' ' << sigma[i] << '\n';
	}
	file_temperature.close();
	cout << "Numerical solution temperature profile in radiative equilibrium stored in file " << fn_temperature_numerical << endl;
	file_irradiance.close();
	cout << "Numerical solution irradiances in radiative equilibrium stored in file " << fn_irradiance_numerical << endl;
	file_errors.close();
	cout << "Absolute errors between normalised numerical and analytical solutions stored in file " << fn_errors << endl;
	
	
	
	/* Stability analysis of numerical solution in radiative equilibrium */
	
	// Prepare variables.
	int n;
	double ddelta;
	
	// Prepare output file.
	char fn_stability[] = DIR_DATA "/stability.dat";
	file_errors << scientific << setprecision(N_PRECISION);
	file_errors.open(fn_stability);
	file_errors << "#N T/T_0 E_U/S_t E_D/S_t" << endl;
	file_errors << "#'1' '1' '1' '1'" << endl;
	
	// Integrate up to ground level.
	for (int i = 0; i < N_STABILITY; i++) {
		n = 1 << i;
		ddelta = (global_delta_g - global_delta_TOA) / n;
		Y[0] = 0.5;
		Y[1] = 1.0;
		Y[2] = 0.0;
		integrate_IVP(n, ddelta, Y, rhs, 3, rungekutta4);
		T_tmp = pow(Y[0], 0.25);
		file_errors << n << ' ' << fabs(T_tmp - temperature_norm(global_delta_g)) << ' ' << fabs(Y[1] - irradiance_upward_norm((global_delta_g))) << ' ' << fabs(Y[2] - irradiance_downward_norm((global_delta_g))) << '\n';
	}
	file_errors.close();
	cout << "Errors for stability analysis of numerical solution stored in file " << fn_stability << endl;
	
	// Tear down.
	delete[] z;
	delete[] P;
	delete[] delta;
	delete[] sigma;
	delete[] T;
	delete[] theta;
	delete[] E_U;
	delete[] E_D;
	
	return 0;
}

void rhs(double t, double const * Y_0, double * R) {
	R[0] = 0.5 * global_D;
	R[1] = global_D * (Y_0[1] - Y_0[0]);
	R[2] = global_D * (Y_0[0] - Y_0[2]);
}
