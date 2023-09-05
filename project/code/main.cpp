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
	double z_TOA; // / m
	double * z, * dz; // / m
	z_TOA = get_altitude(global_P_TOA);
	z = new double[global_N];
	dz = new double[global_N];
	set_layers_uniform(global_z_g, z_TOA, global_N, z, dz);
	
	// Precompute other vertical coordinates.
	double * P; // / Pa
	double * delta;
	P = new double[global_N];
	delta = new double[global_N];
	for (int i = 0; i < global_N; i++) {
		delta[i] = get_optical_depth_z(z[i], global_P_TOA);
		P[i] = get_pressure(z[i]);
	}
	
	// Analytical solution of radiative equilibrium.
	double T_0; // / K
	double * T; // / K
	T_0 = pow((1.0 - global_A) * global_S_0 / 8.0 / global_sigma, 0.25);
	T = new double[global_N];
	ofstream file_plot;
	char filename_plot[] = DIR_DATA "/temperature_radiative_equilibrium.dat";
	file_plot << fixed << setprecision(N_PRECISION);
	file_plot.open(filename_plot);
	file_plot << "#z/z_0 T   delta P   " << endl;
	file_plot << "#''    'K' ''    'Pa'" << endl;
	for (int i = 0; i < global_N; i++) {
		T[i] = pow(1.0 + global_D * delta[i], 0.25);
		file_plot << z[i] / global_z_0 << ' ' << T_0 * T[i] << ' ' << delta[i] << ' ' << P[i] << '\n';
	}
	file_plot.close();
	cout << "Temperature profile calculated, values are stored in file " << filename_plot << endl;
	
	// Clean up.
	delete[] z;
	delete[] dz;
	delete[] P;
	delete[] delta;
	
	return 0;
}
