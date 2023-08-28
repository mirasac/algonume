#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "constants.h"
#include "../../mclib/mclib.h"
#include "radiation.h"
#include "configuration.h" // Load last to redefine some things.

void set_layers_z_uniform(double z_min, double z_max, int n_layers, double z[], double delta_z[]);
void set_absorbers_uniform(int n_layers, int n_absorbers[], absorber_t * absorbers[], absorber_t a1, absorber_t a2);
double get_pressure(double z, double T);
double get_altitude(double P, double T);
double get_sigma(double P, double P_TOA);
double get_theta(double T, double P);

int main(int argc, char * argv[]) {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	
	// Configure atmospheric layers.
	int n_layers = 20;
	double z_TOA; // / m
	double * z, * delta_z; // / m
	z_TOA = 55000.0;
	z = new double[n_layers];
	delta_z = new double[n_layers];
	set_layers_z_uniform(global_z_g, z_TOA, n_layers, z, delta_z);
	
	// Configure absorbers.
	absorber_t H2O, CO2;
	H2O.nu = new double[]{140.0, 1600.0, 3760.0, 5350.0, 7250.0};
	H2O.n_nu = 5;
	H2O.delta_nu = new double[H2O.n_nu]; // MC bandwidths missing.
	CO2.nu = new double[]{667.0, 960.0, 1060.0, 2410.0, 3660.0, 5200.0};
	CO2.n_nu = 6;
	CO2.delta_nu = new double[CO2.n_nu]; // MC bandwidths missing.
	
	int * n_absorbers;
	absorber_t ** absorbers;
	n_absorbers = new int[n_layers];
	absorbers = new absorber_t*[n_layers];
	set_absorbers_uniform(n_layers, n_absorbers, absorbers, H2O, CO2);
	
	// Configure spectral bands.
	int * n_bands;
	double * nu, * delta_nu; // / (1 / cm)
	double dnu, nu_div; // / (1 / m)
	dnu = 10000.0;
	nu_div = spectrum_division_nu();
	// MC continue with three arrays for the bandwidth of each layer.
	
	// Set time integration parameters.
	int n_t = 10680;
	double t, dt, t_min, t_max; // / h
	t_min = 0.0;
	t_max = 85440.0;
	dt = (t_max - t_min) / n_t;
	
	// Initialise temperature output variable.
	double * T; // / K
	T = new double[n_layers];
	for (int i = 0; i < n_layers; i++) {
		T[i] = global_T_earth;
	}
	
	// Prepare support variables.
	double mu; // / rad
	double P, P_TOA; // / Pa
	mu = M_SQRT2 / 2.0;
	
	// Prepare output file.
	ofstream file_plot;
	char filename_plot[] = DIR_DATA "/temperature.dat";
	file_plot.open(filename_plot);
	file_plot << fixed << setprecision(N_PRECISION);
	file_plot << "#t z T P sigma theta" << endl;
	
	// Run model.
	// MC put here initial and boundary conditions.
	P_TOA = get_pressure(z[0], T[0]);
	for (int i_t = 0; i_t <= n_t; i_t++) { // MC switch to while loop to check convergence condition.
		t = t_min + i_t * dt;
		for (int i_z = 0; i_z < n_layers; i_z++) {
			// MC inner integration loop, where radiative calculations and convective adjustment are performed.
			P = get_pressure(z[i_z], T[i_z]);
			file_plot << t << ' ' << z[i_z] << ' ' << T[i_z] << ' ' << P << ' ' << get_sigma(P, P_TOA) << ' ' << get_theta(T[i_z], P) << '\n';
		}
		// MC an additional separate line is needed for values at ground level since they are evaluated separately from the loop on layers.
		file_plot << '\n';
		// MC avoid adding an extra data block at the end of data file.
		if (i_t != n_t) {
			file_plot << '\n';
		}
	}
	file_plot.close();
	//cout << "Temperature profile calculated, values are stored in file " << filename_plot << endl;
	
	// Clean up.
	delete[] z;
	delete[] delta_z;
	delete_absorber(H2O);
	delete_absorber(CO2);
	delete[] n_absorbers;
	delete[] absorbers[0];
	delete[] absorbers;
	
	return 0;
}

void set_layers_z_uniform(double z_min, double z_max, int n_layers, double z[], double delta_z[]) {
	delta_z[0] = (z_max - z_min) / n_layers;
	z[0] = z_max;
	for (int i_z = 1; i_z < n_layers; i_z++) {
		delta_z[i_z] = delta_z[0];
		z[i_z] = z[0] - delta_z[0] * i_z;
	}
}

void set_absorbers_uniform(int n_layers, int n_absorbers[], absorber_t * absorbers[], absorber_t a1, absorber_t a2) {
	n_absorbers[0] = 2;
	absorbers[0] = new absorber_t[n_layers * n_absorbers[0]];
	absorbers[0][0] = a1;
	absorbers[0][1] = a2;
	for (int i = 1; i < n_layers; i++) {
		n_absorbers[i] = n_absorbers[0];
		absorbers[i] = absorbers[0] + i * n_absorbers[0];
		absorbers[i][0] = a1;
		absorbers[i][1] = a2;
	}
}

double get_pressure(double z, double T) {
	return global_P_g * exp(- global_g / global_R_m_air / T * (z - global_z_g));
}

double get_altitude(double P, double T) {
	return global_z_g - global_R_m_air * T / global_g * log(P / global_P_g);
}

double get_sigma(double P, double P_TOA) {
	return (P - P_TOA) / (global_P_g - P_TOA);
}

double get_theta(double T, double P) {
	return T * pow(global_P_0 / P, global_R_m_air / global_c_P_air);
}
