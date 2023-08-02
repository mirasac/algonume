#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "constants.h"
#include "../../mclib/mclib.h"
#include "radiation.h"

void set_layers_z_uniform(double z_TOA, int n_layers, double z[], double dz[]);

int main(int argc, char * argv[]) {
	using namespace std;
	cout << setprecision(N_PRECISION);
	
	// Configure atmospheric layers.
	int const n_layers = 20;
	double z_TOA, z_S; // [m]
	double z[n_layers], dz[n_layers]; // [m]
	z_TOA = 42000.0;
	z_S = 0;
	set_layers_z_uniform(z_TOA, n_layers, z, dz);
	
	return 0;
}

void set_layers_z_uniform(double z_TOA, int n_layers, double z[], double dz[]) {
	dz[0] = z_TOA / n_layers;
	z[0] = z_TOA - dz[0] / 2.0;
	for (int i_z = 1; i_z < n_layers; i_z++) {
		dz[i_z] = dz[0];
		z[i_z] = z[i_z-1] - dz[i_z];
	}
}
