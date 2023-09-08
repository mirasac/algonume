#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "constants.h"
#include "../../mclib/mclib.h"
#include "radiation.h"
#include "utilities.h"
#include "configuration.h" // Include last to allow redefinitions.

int main(int argc, char * argv[]) {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	
	// Configure vertical coordinates.
	double * z, * dz; // / m
	z = new double[global_N];
	dz = new double[global_N];
	set_layers_uniform(global_z_g, get_altitude(global_P_TOA), global_N, z, dz);
	
	// Precompute other vertical coordinates.
	double * P; // / Pa
	double * delta, * sigma;
	delta = new double[global_N];
	P = new double[global_N];
	sigma = new double[global_N];
	for (int i = 0; i < global_N; i++) {
		delta[i] = get_optical_depth_z(z[i], global_P_TOA);
		P[i] = get_pressure(z[i]);
		sigma[i] = get_sigma(P[i], global_P_TOA);
	}
	
	
	
	/* Analytical solution of radiative equilibrium */
	
	double T_0; // / K
	double * T, * theta; // / K
	double S_t, E_U, E_D; // / (W / m^2)
	S_t = (1.0 - global_A) * global_S_0 / 4.0;
	T_0 = pow(S_t / 2.0 / global_sigma, 0.25);
	T = new double[global_N];
	theta = new double[global_N];
	ofstream file_temperature, file_irradiance;
	char fn_temperature_analytical[] = DIR_DATA "/temperature_analytical.dat";
	char fn_irradiance_analytical[] = DIR_DATA "/irradiance_analytical.dat";
	file_temperature << fixed << setprecision(N_PRECISION);
	file_temperature.open(fn_temperature_analytical);
	file_temperature << "#z/z_0  T   delta P    sigma theta" << endl;
	file_temperature << "#'1'    'K' '1'   'Pa' '1'   'K'" << endl;
	file_irradiance << fixed << setprecision(N_PRECISION);
	file_irradiance.open(fn_irradiance_analytical);
	file_irradiance << "#z/z_0  E_U       E_D       delta P    sigma" << endl;
	file_irradiance << "#'1'    'W / m^2' 'W / m^2' '1'   'Pa' '1'" << endl;
	for (int i = 0; i < global_N; i++) {
		T[i] = pow(1.0 + global_D * delta[i], 0.25);
		theta[i] = get_theta(T[i], P[i]);
		file_temperature << z[i] / global_z_0 << ' ' << T_0 * T[i] << ' ' << delta[i] << ' ' << P[i] << ' ' << sigma[i] << ' ' << theta[i] << '\n';
		E_U = 0.5 * (2.0 + global_D * delta[i]);
		E_D = 0.5 * global_D * delta[i];
		file_irradiance << z[i] / global_z_0 << ' ' << E_U << ' ' << E_D << ' ' << delta[i] << ' ' << P[i] << ' ' << sigma[i] << '\n';
	}
	file_temperature.close();
	cout << "Analytical solution temperature profile in radiative equilibrium calculated, values are stored in file " << fn_temperature_analytical << endl;
	file_irradiance.close();
	cout << "Analytical solution irradiances in radiative equilibrium calculated, values are stored in file " << fn_irradiance_analytical << endl;
	
	// Tear down.
	delete[] z;
	delete[] dz;
	delete[] P;
	delete[] delta;
	delete[] sigma;
	delete[] T;
	delete[] theta;
	
	return 0;
}
