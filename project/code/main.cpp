#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "constants.h"
#include "../../mclib/mclib.h"
#include "radiation.h"

void set_layers_z_uniform(double z_S, double z_TOA, int n_layers, double z[], double delta_z[]);
void set_absorbers_default(int n_layers, int n_absorbers[], absorber_t * absorbers[], absorber_t a1, absorber_t a2);

int main(int argc, char * argv[]) {
	using namespace std;
	cout << setprecision(N_PRECISION);
	
	// Configure atmospheric layers.
	int n_layers = 20;
	double z_TOA, z_S; // [m]
	double * z, * delta_z; // [m]
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
	
	// MC debug Planck function.
	ofstream plot_file;
	plot_file.open("MC.dat");
	plot_file << fixed << setprecision(N_PRECISION);
	plot_file << "#nu B(nu)" << endl;
	for (int i = 0; i < CO2.n_nu; i++) {
		plot_file << CO2.nu[i] << ' ' << planck_law(100.0 * CO2.nu[i], global_T_0) << endl;
	}
	plot_file.close();
	
	int * n_absorbers;
	absorber_t ** absorbers;
	n_absorbers = new int[n_layers];
	absorbers = new absorber_t*[n_layers];
	set_absorbers_default(n_layers, n_absorbers, absorbers, H2O, CO2);
	
	// Configure spectral bands.
	int * n_bands;
	double * nu, * delta_nu; // [1 / cm]
	double dnu = 100.0; // [1 / cm]
	// MC continue with three arrays for the bandwidth of each layer.
	
	// Clean up.
	delete[] z;
	delete[] delta_z;
	delete_absorber(H2O);
	delete_absorber(CO2);
	delete[] n_absorbers;
	for (int i = 0; i < n_layers; i++) {
		delete[] absorbers[i];
	}
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

void set_absorbers_default(int n_layers, int n_absorbers[], absorber_t * absorbers[], absorber_t a1, absorber_t a2) {
	for (int i_layers = 0; i_layers < n_layers; i_layers++) {
		absorbers[i_layers] = new absorber_t[]{a1, a2};
		n_absorbers[i_layers] = 2;
	}
}
