#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "constants.h"
#include "convection.h"
#include "../../mclib/mclib.h"
#include "radiation.h"
#include "utilities.h"
#include "configuration.h" // Include last to allow redefinitions.

double * global_z; // / m
double * global_P; // / Pa
double * global_delta, * global_sigma;
double global_T_0; // / K
double global_Delta_t; // / s

void rhs_delta(double t, double const * Y_0, double * R);

void rhs_t(double t, double const * Y_0, double * R);

int main(int argc, char * argv[]) {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	
	// Configure vertical coordinates.
	double z_TOA, dz; // / m
	z_TOA = get_altitude(global_P_TOA);
	dz = (z_TOA - const_z_g) / global_N;
	global_z = new double[global_N + 1];
	global_delta = new double[global_N + 1];
	global_P = new double[global_N + 1];
	global_sigma = new double[global_N + 1];
	for (int i = 0; i <= global_N; i++) {
		global_z[i] = z_TOA - i * dz;
		global_delta[i] = get_optical_depth_z(global_z[i], global_P_TOA);
		global_P[i] = get_pressure(global_z[i]);
		global_sigma[i] = get_sigma(global_P[i], global_P_TOA);
	}
	
	
	
	/* Analytical solution in radiative equilibrium */
	
	cout << "Analytical solution in radiative equilibrium" << endl;
	
	// Prepare variables.
	double * Y_3_ana, * theta_norm, * Y_1_ana, * Y_2_ana;
	Y_3_ana = new double[global_N + 1];
	theta_norm = new double[global_N + 1];
	Y_1_ana = new double[global_N + 1];
	Y_2_ana = new double[global_N + 1];
	global_T_0 = pow(global_S_t / const_sigma, 0.25);
	
	// Prepare output files.
	ofstream file_temperature, file_irradiance;
	char fn_temperature_analytical[] = DIR_DATA "/temperature_analytical.dat";
	char fn_irradiance_analytical[] = DIR_DATA "/irradiance_analytical.dat";
	file_temperature << fixed << setprecision(N_PRECISION);
	file_temperature.open(fn_temperature_analytical);
	file_temperature << "#z/z_0 T/T_0 delta P sigma theta/T_0" << endl;
	file_temperature << "#'1' '1' '1' 'Pa' '1' '1'" << endl;
	file_irradiance << fixed << setprecision(N_PRECISION);
	file_irradiance.open(fn_irradiance_analytical);
	file_irradiance << "#z/z_0 E_U/S_t E_D/S_t delta P sigma" << endl;
	file_irradiance << "#'1' '1' '1' '1' 'Pa' '1'" << endl;
	
	// Plot analytical solutions.
	for (int i = 0; i <= global_N; i++) {
		Y_3_ana[i] = temperature_norm(global_delta[i]);
		theta_norm[i] = get_theta(Y_3_ana[i], global_P[i]);
		file_temperature << global_z[i] / global_z_0 << ' '
			<< Y_3_ana[i] << ' '
			<< global_delta[i] << ' '
			<< global_P[i] << ' '
			<< global_sigma[i] << ' '
			<< theta_norm[i] << '\n';
		Y_1_ana[i] = irradiance_upward_norm(global_delta[i]);
		Y_2_ana[i] = irradiance_downward_norm(global_delta[i]);
		file_irradiance << global_z[i] / global_z_0 << ' '
			<< Y_1_ana[i] << ' '
			<< Y_2_ana[i] << ' '
			<< global_delta[i] << ' '
			<< global_P[i] << ' '
			<< global_sigma[i] << '\n';
	}
	file_temperature.close();
	cout << "- Temperature stored in file " << fn_temperature_analytical << endl;
	file_irradiance.close();
	cout << "- Irradiances stored in file " << fn_irradiance_analytical << endl;
	
	
	
	/* Numerical solution in radiative equilibrium */
	
	cout << "Numerical solution in radiative equilibrium" << endl;
	
	// Prepare variables.
	int i_t, is_steady;
	double * Y_3, * Y_3_prev, * Y_1, * Y_2;
	double t; // / s
	double Y[3];
	double Y_3_tmp, theta_tmp, theta_PDE_tmp;
	Y_3 = new double[3 * (global_N + 1)]; // Use to store irradiances.
	Y_3_prev = new double[global_N + 1];
	Y_1 = Y_3 + global_N + 1;
	Y_2 = Y_1 + global_N + 1;
	global_Delta_t = 10 * 24 * 3600.0;
	
	// Set initial values PDEs.
	i_t = 0;
	Y_3[0] = temperature_norm(const_delta_TOA);
	for (int i = 1; i <= global_N; i++) {
		Y_3[i] = const_T_g / global_T_0;
	}
	Y_1[0] = 1.0;
	Y_2[0] = 0.0;
	
	// Integrate PDEs.
	do {
		Y_3_prev[0] = Y_3[0];
		i_t++;
		t = i_t * global_Delta_t;
		Y[1] = 1.0;
		Y[2] = 0.0;
		for (int i = 1; i <= global_N; i++) {
			Y_3_prev[i] = Y_3[i];
			Y[0] = Y_3[i]*Y_3[i]*Y_3[i]*Y_3[i];
			rungekutta4(global_delta[i], global_delta[i] - global_delta[i-1], Y, rhs_delta, 3);
			Y_1[i] = Y[1];
			Y_2[i] = Y[2];
		}
		eulerstep(t, global_Delta_t, Y_3, rhs_t, global_N + 1);
		// Temperature profile steady state condition.
		is_steady = 1;
		for (int i = 0; i <= global_N; i++) {
			if (fabs(Y_3[i] - Y_3_prev[i]) > TOLERANCE) {
				is_steady = 0;
				break;
			}
		}
	} while (! is_steady);
	cout << "- Steady state reached in " << i_t << " iterations" << endl;
	
	// Prepare output files.
	char fn_temperature_numerical[] = DIR_DATA "/temperature_numerical.dat";
	char fn_irradiance_numerical[] = DIR_DATA "/irradiance_numerical.dat";
	file_temperature << fixed << setprecision(N_PRECISION);
	file_temperature.open(fn_temperature_numerical);
	file_temperature << "#z/z_0 T/T_0 T_PDE/t_0 T_err/T_0 T_PDE_err/T_0 delta P sigma theta/T_0 theta_PDE/T_0 theta_err/T_0 theta_PDE_err/T_0" << endl;
	file_temperature << "#'1' '1' '1' '1' '1' '1' 'Pa' '1' '1' '1' '1' '1'" << endl;
	file_irradiance << fixed << setprecision(N_PRECISION);
	file_irradiance.open(fn_irradiance_numerical);
	file_irradiance << "#z/z_0 E_U/S_t E_U_PDE/S_t E_D/S_t E_D_PDE/S_t E_U_err/S_t E_U_PDE_err/S_t E_D_err/S_t E_D_PDE_err/S_t delta P sigma" << endl;
	file_irradiance << "#'1' '1' '1' '1' '1' '1' '1' '1' '1' '1' 'Pa' '1'" << endl;
	
	// Set initial values ODEs.
	Y[0] = 0.5;
	Y[1] = 1.0;
	Y[2] = 0.0;
	Y_3_tmp = pow(Y[0], 0.25);
	theta_tmp = get_theta(Y_3_tmp, global_P[0]);
	theta_PDE_tmp = get_theta(Y_3[0], global_P[0]);
	file_temperature << global_z[0] / global_z_0 << ' '
		<< Y_3_tmp << ' '
		<< Y_3[0] << ' '
		<< scientific << fabs(Y_3_tmp - Y_3_ana[0]) << ' '
		<< fabs(Y_3[0] - Y_3_ana[0]) << fixed << ' '
		<< global_delta[0] << ' '
		<< global_P[0] << ' '
		<< global_sigma[0] << ' '
		<< theta_tmp << ' '
		<< theta_PDE_tmp << ' '
		<< scientific << fabs(theta_tmp - theta_norm[0]) << ' '
		<< fabs(theta_PDE_tmp - theta_norm[0]) << fixed << '\n';
	file_irradiance << global_z[0] / global_z_0 << ' '
		<< Y[1] << ' '
		<< Y_1[0] << ' '
		<< Y[2] << ' '
		<< Y_2[0] << ' '
		<< scientific << fabs(Y[1] - Y_1_ana[0]) << ' '
		<< fabs(Y_1[0] - Y_1_ana[0]) << ' '
		<< fabs(Y[2] - Y_2_ana[0]) << ' '
		<< fabs(Y_2[0] - Y_2_ana[0]) << fixed << ' '
		<< global_delta[0] << ' '
		<< global_P[0] << ' '
		<< global_sigma[0] << '\n';
	
	// Integrate ODEs.
	for (int i = 1; i <= global_N; i++) {
		rungekutta4(global_delta[i], global_delta[i] - global_delta[i-1], Y, rhs_delta, 3);
		Y_3_tmp = pow(Y[0], 0.25);
		theta_tmp = get_theta(Y_3_tmp, global_P[i]);
		theta_PDE_tmp = get_theta(Y_3[i], global_P[i]);
		file_temperature << global_z[i] / global_z_0 << ' '
			<< Y_3_tmp << ' '
			<< Y_3[i] << ' '
			<< scientific << fabs(Y_3_tmp - Y_3_ana[i]) << ' '
			<< fabs(Y_3[i] - Y_3_ana[i]) << fixed << ' '
			<< global_delta[i] << ' '
			<< global_P[i] << ' '
			<< global_sigma[i] << ' '
			<< theta_tmp << ' '
			<< theta_PDE_tmp << ' '
			<< scientific << fabs(theta_tmp - theta_norm[i]) << ' '
			<< fabs(theta_PDE_tmp - theta_norm[i]) << fixed << '\n';
		file_irradiance << global_z[i] / global_z_0 << ' '
			<< Y[1] << ' '
			<< Y_1[i] << ' '
			<< Y[2] << ' '
			<< Y_2[i] << ' '
			<< scientific << fabs(Y[1] - Y_1_ana[i]) << ' '
			<< fabs(Y_1[i] - Y_1_ana[i]) << ' '
			<< fabs(Y[2] - Y_2_ana[i]) << ' '
			<< fabs(Y_2[i] - Y_2_ana[i]) << fixed << ' '
			<< global_delta[i] << ' '
			<< global_P[i] << ' '
			<< global_sigma[i] << '\n';
	}
	file_temperature.close();
	cout << "- Temperature stored in file " << fn_temperature_numerical << endl;
	file_irradiance.close();
	cout << "- Irradiances stored in file " << fn_irradiance_numerical << endl;
	
	
	
	/* Stability analysis of numerical solution in radiative equilibrium */
	
	cout << "Stability analysis of numerical solution in radiative equilibrium" << endl;
	
	// Prepare output file.
	ofstream file_stability;
	char fn_stability[] = DIR_DATA "/stability.dat";
	file_stability << scientific << setprecision(N_PRECISION);
	file_stability.open(fn_stability);
	file_stability << "#N T_err(delta_g)/T_0 E_U_err(delta_g)/S_t E_D_err(delta_g)/S_t" << endl;
	file_stability << "#'1' '1' '1' '1'" << endl;
	
	// Integrate up to ground level.
	int n;
	for (int i = 0; i < N_STABILITY; i++) {
		n = 1 << i;
		Y[0] = 0.5;
		Y[1] = 1.0;
		Y[2] = 0.0;
		integrate_IVP(n, (global_delta_g - const_delta_TOA) / n, Y, rhs_delta, 3, rungekutta4);
		Y_3_tmp = pow(Y[0], 0.25);
		file_stability << n << ' '
			<< fabs(Y_3_tmp - Y_3_ana[global_N]) << ' '
			<< fabs(Y[1] - Y_1_ana[global_N]) << ' '
			<< fabs(Y[2] - Y_2_ana[global_N]) << '\n';
	}
	file_stability.close();
	cout << "- Errors stored in file " << fn_stability << endl;
	
	
	
	/* Radiative-convective equilibrium */
	
	cout << "Radiative-convective equilibrium" << endl;
	
	// Set initial values.
	i_t = 0;
	Y_3[0] = temperature_norm(const_delta_TOA);
	for (int i = 1; i <= global_N; i++) {
		Y_3[i] = const_T_g / global_T_0;
	}
	Y_1[0] = 1.0;
	Y_2[0] = 0.0;
	
	// Run model.
	do {
		Y_3_prev[0] = Y_3[0];
		i_t++;
		t = i_t * global_Delta_t;
		Y[1] = 1.0;
		Y[2] = 0.0;
		for (int i = 1; i <= global_N; i++) {
			Y_3_prev[i] = Y_3[i];
			Y[0] = Y_3[i]*Y_3[i]*Y_3[i]*Y_3[i];
			rungekutta4(global_delta[i], global_delta[i] - global_delta[i-1], Y, rhs_delta, 3);
			Y_1[i] = Y[1];
			Y_2[i] = Y[2];
		}
		eulerstep(t, global_Delta_t, Y_3, rhs_t, global_N + 1);
		convective_adjustment(const_Gamma_0 / global_T_0, global_N, Y_3, global_z);
		// Temperature profile steady state condition.
		is_steady = 1;
		for (int i = 0; i <= global_N; i++) {
			if (fabs(Y_3[i] - Y_3_prev[i]) > TOLERANCE) {
				is_steady = 0;
				break;
			}
		}
	} while (! is_steady);
	cout << "- Steady state reached in " << i_t << " iterations" << endl;
	
	// Prepare output files.
	char fn_temperature_RCM[] = DIR_DATA "/temperature_RCM.dat";
	char fn_irradiance_RCM[] = DIR_DATA "/irradiance_RCM.dat";
	file_temperature << fixed << setprecision(N_PRECISION);
	file_temperature.open(fn_temperature_RCM);
	file_temperature << "#z/z_0 T/T_0 T_err/T_0 delta P sigma theta/T_0 theta_err/T_0" << endl;
	file_temperature << "#'1' '1' '1' '1' 'Pa' '1' '1' '1'" << endl;
	file_irradiance << fixed << setprecision(N_PRECISION);
	file_irradiance.open(fn_irradiance_RCM);
	file_irradiance << "#z/z_0 E_U/S_t E_D/S_t E_U_err/S_t E_D_err/S_t delta P sigma" << endl;
	file_irradiance << "#'1' '1' '1' '1' '1' '1' 'Pa' '1'" << endl;
	
	// Print output values.
	for (int i = 0; i <= global_N; i++) {
		theta_tmp = get_theta(Y_3[i], global_P[i]);
		file_temperature << global_z[i] / global_z_0 << ' '
			<< Y_3[i] << ' '
			<< scientific << fabs(Y_3[i] - Y_3_ana[i]) << fixed << ' '
			<< global_delta[i] << ' '
			<< global_P[i] << ' '
			<< global_sigma[i] << ' '
			<< theta_tmp << ' '
			<< scientific << fabs(theta_tmp - theta_norm[i]) << fixed << '\n';
		file_irradiance << global_z[i] / global_z_0 << ' '
			<< Y_1[i] << ' '
			<< Y_2[i] << ' '
			<< scientific << fabs(Y_1[i] - Y_1_ana[i]) << ' '
			<< fabs(Y_2[i] - Y_2_ana[i]) << fixed << ' '
			<< global_delta[i] << ' '
			<< global_P[i] << ' '
			<< global_sigma[i] << '\n';
	}
	file_temperature.close();
	cout << "- Temperature stored in file " << fn_temperature_RCM << endl;
	file_irradiance.close();
	cout << "- Irradiances stored in file " << fn_irradiance_RCM << endl;
	
	// Tear down.
	delete[] global_z;
	delete[] global_P;
	delete[] global_delta;
	delete[] global_sigma;
	delete[] Y_3_ana;
	delete[] theta_norm;
	delete[] Y_1_ana;
	delete[] Y_2_ana;
	delete[] Y_3;
	delete[] Y_3_prev;
	
	return 0;
}

void rhs_delta(double t, double const * Y_0, double * R) {
	R[0] = 0.5 * const_D;
	R[1] = const_D * (Y_0[1] - Y_0[0]);
	R[2] = const_D * (Y_0[0] - Y_0[2]);
}

void rhs_t(double t, double const * Y_0, double * R) {
	double const * Y_1, * Y_2;
	Y_1 = Y_0 + global_N + 1;
	Y_2 = Y_1 + global_N + 1;
	R[0] = 0.0;
	for (int i = 1; i <= global_N; i++) {
		R[i] = global_S_t * global_mu_m / (const_c_P * global_T_0) * (Y_1[i] - Y_1[i-1] - Y_2[i] + Y_2[i-1]) / (global_delta[i] - global_delta[i-1]);
	}
}
