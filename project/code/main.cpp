#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "constants.h"
#include "../../mclib/mclib.h"
#include "radiation.h"
#include "utilities.h"
#include "configuration.h" // Include last to allow redefinitions.

void rhs_Y(double t, double const * Y_0, double * R);

void rhs_T(double t, double const * Y_0, double * R);

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
	double * T_norm, * theta_norm, * Y_1, * Y_2;
	double T_0; // / K
	T_norm = new double[global_N + 1];
	theta_norm = new double[global_N + 1];
	Y_1 = new double[global_N + 1];
	Y_2 = new double[global_N + 1];
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
		T_norm[i] = temperature_norm(delta[i]);
		theta_norm[i] = get_theta(T_norm[i], P[i]);
		file_temperature << z[i] / global_z_0 << ' ' << T_0 * T_norm[i] << ' ' << delta[i] << ' ' << P[i] << ' ' << sigma[i] << ' ' << T_0 * theta_norm[i] << '\n';
		Y_1[i] = irradiance_upward_norm(delta[i]);
		Y_2[i] = irradiance_downward_norm(delta[i]);
		file_irradiance << z[i] / global_z_0 << ' ' << Y_1[i] << ' ' << Y_2[i] << ' ' << delta[i] << ' ' << P[i] << ' ' << sigma[i] << '\n';
	}
	file_temperature.close();
	cout << "Analytical solution temperature profile in radiative equilibrium calculated, values are stored in file " << fn_temperature_analytical << endl;
	file_irradiance.close();
	cout << "Analytical solution irradiances in radiative equilibrium calculated, values are stored in file " << fn_irradiance_analytical << endl;
	
	
	
	/* Numerical solution in radiative equilibrium */
	
	// Prepare variables and set initial values.
	double Y[3];
	double T_tmp, theta_tmp, ddelta;
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
	file_temperature << z[0] / global_z_0 << ' ' << T_0 * T_tmp << ' ' << delta[0] << ' ' << P[0] << ' ' << sigma[0] << ' ' << T_0 * theta_norm[0] << '\n';
	file_irradiance << z[0] / global_z_0 << ' ' << Y[1] << ' ' << Y[2] << ' ' << delta[0] << ' ' << P[0] << ' ' << sigma[0] << '\n';
	file_errors << z[0] / global_z_0 << ' ' << fabs(T_tmp - T_norm[0]) << ' ' << fabs(theta_tmp - theta_norm[0]) << ' ' << fabs(Y[1] - Y_1[0]) << ' ' << fabs(Y[2] - Y_2[0]) << ' ' << delta[0] << ' ' << P[0] << ' ' << sigma[0] << '\n';
	for (int i = 1; i <= global_N; i++) {
		ddelta = delta[i] - delta[i-1];
		rungekutta4(delta[i], ddelta, Y, rhs_Y, 3);
		T_tmp = pow(Y[0], 0.25);
		theta_tmp = get_theta(T_tmp, P[i]);
		file_temperature << z[i] / global_z_0 << ' ' << T_0 * T_tmp << ' ' << delta[i] << ' ' << P[i] << ' ' << sigma[i] << ' ' << T_0 * theta_norm[i] << '\n';
		file_irradiance << z[i] / global_z_0 << ' ' << Y[1] << ' ' << Y[2] << ' ' << delta[i] << ' ' << P[i] << ' ' << sigma[i] << '\n';
		file_errors << z[i] / global_z_0 << ' ' << fabs(T_tmp - T_norm[i]) << ' ' << fabs(theta_tmp - theta_norm[i]) << ' ' << fabs(Y[1] - Y_1[i]) << ' ' << fabs(Y[2] - Y_2[i]) << ' ' << delta[i] << ' ' << P[i] << ' ' << sigma[i] << '\n';
	}
	file_temperature.close();
	cout << "Numerical solution temperature profile in radiative equilibrium stored in file " << fn_temperature_numerical << endl;
	file_irradiance.close();
	cout << "Numerical solution irradiances in radiative equilibrium stored in file " << fn_irradiance_numerical << endl;
	file_errors.close();
	cout << "Absolute errors between normalised numerical and analytical solutions stored in file " << fn_errors << endl;
	
	
	
	/* Stability analysis of numerical solution in radiative equilibrium */
	
	// Prepare output file.
	char fn_stability[] = DIR_DATA "/stability.dat";
	file_errors << scientific << setprecision(N_PRECISION);
	file_errors.open(fn_stability);
	file_errors << "#N T/T_0 E_U/S_t E_D/S_t" << endl;
	file_errors << "#'1' '1' '1' '1'" << endl;
	
	// Integrate up to ground level.
	int n;
	for (int i = 0; i < N_STABILITY; i++) {
		n = 1 << i;
		ddelta = (global_delta_g - global_delta_TOA) / n;
		Y[0] = 0.5;
		Y[1] = 1.0;
		Y[2] = 0.0;
		integrate_IVP(n, ddelta, Y, rhs_Y, 3, rungekutta4);
		T_tmp = pow(Y[0], 0.25);
		file_errors << n << ' ' << fabs(T_tmp - temperature_norm(global_delta_g)) << ' ' << fabs(Y[1] - irradiance_upward_norm((global_delta_g))) << ' ' << fabs(Y[2] - irradiance_downward_norm((global_delta_g))) << '\n';
	}
	file_errors.close();
	cout << "Errors for stability analysis of numerical solution stored in file " << fn_stability << endl;
	
	
	
	/* Numerical solution in radiative-convective equilibrium */
	
	// Set integration parameters.
	int i_t, n_T, i_0_Y_1, i_0_Y_2, i_0_ddelta, is_steady;
	double * T, * T_prev;
	double t, dt; // / s
	i_t = 0;
	i_0_Y_1 = global_N + 1;
	i_0_Y_2 = i_0_Y_1 + global_N + 1;
	i_0_ddelta = i_0_Y_2 + global_N + 1;
	n_T = i_0_ddelta + global_N + 1;
	T = new double[n_T]; // Used also to store irradiances and ddelta.
	T_prev = new double[global_N + 1]; // Used also to store irradiances and ddelta.
	dt = 12 * 3600.0;
	
	// Set initial values.
	T[0] = T_0 * temperature_norm(global_delta_TOA);
	for (int i = 1; i <= global_N; i++) {
		T[i] = global_T_g;
	}
	T[i_0_Y_1] = 1.0;
	T[i_0_Y_2] = 0.0;
	T[i_0_ddelta] = 0.0;
	
	// Run model.
	do {
		T_prev[0] = T[0];
		t = i_t * dt;
		Y[1] = 1.0;
		Y[2] = 0.0;
		for (int i = 1; i <= global_N; i++) {
			T_prev[i] = T[i];
			ddelta = delta[i] - delta[i-1];
			Y[0] = global_sigma * T[i]*T[i]*T[i]*T[i] / global_S_t;
			rungekutta4(delta[i], ddelta, Y, rhs_Y, 3);
			T[i_0_Y_1 + i] = Y[1];
			T[i_0_Y_2 + i] = Y[2];
			T[i_0_ddelta + i] = ddelta;
		}
		eulerstep(t, dt, T, rhs_T, global_N + 1);
		// MC add here convective adjustment.
		is_steady = 1;
		for (int i = 0; i <= global_N;  i++) {
			if (fabs(T[i] - T_prev[i]) > TOLERANCE) {
				is_steady = 0;
				break;
			}
		}
		i_t++;
	} while (! is_steady); // Temperature profile steady state.
	cout << "Steady state reached in " << i_t << " iterations" << endl;
	
	// Prepare output files.
	char fn_temperature_steady[] = DIR_DATA "/temperature_steady.dat";
	char fn_irradiance_steady[] = DIR_DATA "/irradiance_steady.dat";
	char fn_errors_steady[] = DIR_DATA "/errors_steady.dat";
	file_temperature << fixed << setprecision(N_PRECISION);
	file_temperature.open(fn_temperature_steady);
	file_temperature << "#z/z_0 T T_err/T_0 delta P sigma theta theta_err/T_0" << endl;
	file_temperature << "#'1' 'K' '1' '1' 'Pa' '1' 'K' '1'" << endl;
	file_irradiance << fixed << setprecision(N_PRECISION);
	file_irradiance.open(fn_irradiance_steady);
	file_irradiance << "#z/z_0 E_U/S_t E_D/S_t E_U_err/S_t E_D_err/S_t delta P sigma" << endl;
	file_irradiance << "#'1' '1' '1' '1' '1' '1' 'Pa' '1'" << endl;
	file_errors << scientific << setprecision(N_PRECISION);
	file_errors.open(fn_errors_steady);
	file_errors << "#z/z_0 T theta E_U/S_t E_D/S_t delta P sigma" << endl;
	file_errors << "#'1' 'K' 'K' '1' '1' '1' 'Pa' '1'" << endl;
	
	// Print output values.
	for (int i = 1; i <= global_N; i++) {
		theta_tmp = get_theta(T[i], P[i]);
		file_temperature << z[i] / global_z_0 << ' ' << T[i] << ' ' << delta[i] << ' ' << P[i] << ' ' << sigma[i] << ' ' << theta_tmp << '\n';
		file_irradiance << z[i] / global_z_0 << ' ' << T[i_0_Y_1 + i] << ' ' << T[i_0_Y_2 + i] << ' ' << delta[i] << ' ' << P[i] << ' ' << sigma[i] << '\n';
		file_errors << z[i] / global_z_0 << ' ' << fabs(T[i] - T_0 * T_norm[i]) << ' ' << fabs(theta_tmp - T_0 * theta_norm[i]) << ' ' << fabs(T[i_0_Y_1 + i] - Y_1[i]) << ' ' << fabs(T[i_0_Y_2 + i] - Y_2[i]) << ' ' << delta[i] << ' ' << P[i] << ' ' << sigma[i] << '\n';
	}
	file_temperature.close();
	cout << "Numerical solution temperature profile in radiative-convective equilibrium stored in file " << fn_temperature_steady << endl;
	file_irradiance.close();
	cout << "Numerical solution irradiances in radiative-convective equilibrium stored in file " << fn_irradiance_steady << endl;
	file_errors.close();
	cout << "Absolute errors between normalised numerical and analytical solutions from the PDE integration stored in file " << fn_errors_steady << endl;
	
	// Tear down.
	delete[] z;
	delete[] P;
	delete[] delta;
	delete[] sigma;
	delete[] T_norm;
	delete[] theta_norm;
	delete[] Y_1;
	delete[] Y_2;
	delete[] T;
	delete[] T_prev;
	
	return 0;
}

void rhs_Y(double t, double const * Y_0, double * R) {
	R[0] = 0.5 * global_D;
	R[1] = global_D * (Y_0[1] - Y_0[0]);
	R[2] = global_D * (Y_0[0] - Y_0[2]);
}

void rhs_T(double t, double const * Y_0, double * R) {
	double const * Y_1, * Y_2, * ddelta;
	Y_1 = Y_0 + global_N + 1;
	Y_2 = Y_1 + global_N + 1;
	ddelta = Y_2 + global_N + 1;
	R[0] = 0.0;
	for (int i = 1; i <= global_N; i++) {
		R[i] = global_S_t * global_mu_m * global_D / global_c_P * (Y_1[i] - Y_1[i-1] - Y_2[i] + Y_2[i-1]) / ddelta[i];
	}
}
