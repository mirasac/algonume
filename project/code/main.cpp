#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "constants.h"
#include "../../mclib/mclib.h"
#include "radiation.h"

void set_layers_z_uniform(double z_TOA, int n_layers, double z[], double delta_z[]);

int main(int argc, char * argv[]) {
	using namespace std;
	cout << setprecision(N_PRECISION);
	
	// Configure atmospheric layers.
	int n_layers = 20;
	double z_TOA, z_S; // [m]
	double * z, * delta_z; // [m]
	z_TOA = 42000.0;
	z_S = 0;
	z = new double[n_layers];
	delta_z = new double[n_layers];
	set_layers_z_uniform(z_TOA, n_layers, z, delta_z);
	
	// Configure absorbers.
	absorber_t O3, CO2, soot;
	// MC continue with detailing the optical bandwidths of absorbers.
	
	int n_absorbers[n_layers];
	absorber_t ** absorbers;
	absorbers = new absorber_t*[n_layers];
	for (int i_layers = 0; i_layers < n_layers; i_layers++) {
		cout << z[i_layers] << endl; // MC debug;
		n_absorbers[i_layers] = 2;
		absorbers[i_layers] = new absorber_t[n_absorbers[i_layers]];
		for (int i_a = 0; i_a < n_absorbers[i_layers]; i_a++) {
			absorbers[i_layers][i_a] = CO2;
		}
	}
	
	// Clean up.
	delete[] z;
	delete[] delta_z;
	for (int i = 0; i < n_layers; i++) {
		delete[] absorbers[i];
	}
	delete[] absorbers;
	
	return 0;
}

void set_layers_z_uniform(double z_TOA, int n_layers, double z[], double delta_z[]) {
	delta_z[0] = z_TOA / n_layers;
	z[0] = z_TOA - delta_z[0] / 2.0;
	for (int i_z = 1; i_z < n_layers; i_z++) {
		delta_z[i_z] = delta_z[0];
		z[i_z] = z[0] - delta_z[0] * i_z;
	}
}
