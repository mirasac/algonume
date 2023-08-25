#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "constants.h"
#include "../../mclib/mclib.h"
#include "radiation.h"
#include "configuration.h" // Load last to redefine some things.

void set_layers_z_uniform(double z_S, double z_TOA, int n_layers, double z[], double delta_z[]);
void set_absorbers_uniform(int n_layers, int n_absorbers[], absorber_t * absorbers[], absorber_t a1, absorber_t a2);

int main(int argc, char * argv[]) {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	
	// Configure atmospheric layers.
	int n_layers = 20;
	double z_TOA, z_S; // / m
	double * z, * delta_z; // / m
	z_TOA = 42000.0;
	z_S = 0.0;
	z = new double[n_layers];
	delta_z = new double[n_layers];
	set_layers_z_uniform(z_S, z_TOA, n_layers, z, delta_z);
	
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
	
	// MC draft variables, decide where to put them.
	double mu = M_SQRT2 / 2.0; // / rad
	
	// Prepare output file.
	ofstream file_plot;
	char filename_plot[] = DIR_DATA "/temperature.dat";
	file_plot.open(filename_plot);
	file_plot << fixed << setprecision(N_PRECISION);
	file_plot << "#t P T" << endl;
	
	// Run model.
	for (int i_t = 0; i_t <= n_t; i_t++) {
		t = t_min + i_t * dt;
		// MC here put inner integration loop, where radiative calculations and convective adjustment are performed.
		file_plot << '\n' << endl; // MC put double '\n' to reduce write time.
	}
	file_plot.close();
	cout << "Temperature profile calculated, values are stored in file " << filename_plot << endl;
	
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

void set_layers_z_uniform(double z_S, double z_TOA, int n_layers, double z[], double delta_z[]) {
	delta_z[0] = (z_TOA - z_S) / n_layers;
	z[0] = z_TOA - delta_z[0] / 2.0;
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
