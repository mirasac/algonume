#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "constants.h"
#include "../../mclib/mclib.h"
#include "radiation.h"
#include "utilities.h"
#include "configuration.h" // Include last to allow redefinitions.

void rhs(double t, const double * Y_0, double * R, int const n_eq);

int main(int argc, char * argv[]) {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	
	// Configure atmospheric layers.
	int n_layers;
	double * z, * delta_z; // / m
	n_layers = 20;
	z = new double[n_layers];
	delta_z = new double[n_layers];
	set_layers_z_uniform(global_z_g, global_z_TOA, n_layers, z, delta_z);
	
	// Configure absorbers.
	absorber_t H2O, CO2;
	H2O.nu = new double[]{14000.0, 160000.0, 376000.0, 535000.0, 725000.0};
	H2O.n_nu = 5;
	H2O.delta_nu = new double[H2O.n_nu]; // MC bandwidths missing.
	CO2.nu = new double[]{66700.0, 96000.0, 106000.0, 241000.0, 366000.0, 520000.0};
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
	int i_t;
	double t, dt; // / s
	dt = 8 * 3600.0;
	
	// Initialise temperature output variable, set initial and boundary conditions.
	double T_TOA; // / K
	double * T; // / K
	T = new double[2 * n_layers + 3]; // Use also for parameters.
	for (int i = 0; i < n_layers; i++) {
		T[i] = global_T_earth;
		T[n_layers + 1 + i] = z[i];
	}
	T[n_layers] = global_T_earth; // Temperature at ground level.
	T[2 * n_layers + 1] = global_z_g;
	T[2 * n_layers + 2] = nu_div;
	
	// Prepare support variables.
	double mu;
	double P, P_TOA; // / Pa
	mu = M_SQRT2 / 2.0;
	
	// Prepare output file.
	ofstream file_plot;
	char filename_plot[] = DIR_DATA "/temperature.dat";
	file_plot.open(filename_plot);
	file_plot << fixed << setprecision(N_PRECISION);
	file_plot << "#t   z   T   P    sigma theta" << endl;
	file_plot << "#(h) (m) (K) (Pa) ()    (K)";
	
	// Run model.
	i_t = 0;
	do {
		file_plot << '\n';
		t = i_t * dt;
		T_TOA = T[0];
		P_TOA = get_pressure(z[0], T_TOA);
		for (int i_z = 0; i_z < n_layers; i_z++) {
			file_plot << t << ' ' << z[i_z] << ' ' << T[i_z] << ' ' << P << ' ' << get_sigma(P, P_TOA) << ' ' << get_theta(T[i_z], P) << '\n';
		}
		file_plot << t << ' ' << global_z_g << ' ' << T[n_layers] << ' ' << P << ' ' << get_sigma(P, P_TOA) << ' ' << get_theta(T[n_layers], P) << '\n';
		eulerstep(t, dt, T, rhs, n_layers);
		// MC an additional separate line of calculation is needed for values at ground level.
		file_plot << '\n';
		i_t++;
	} while (fabs(T[0] - T_TOA) > TOLERANCE); // Equilibrium condition at TOA.
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

void rhs(double t, double const * Y_0, double * R, int const n_eq) {
	int i_0_z, i_0_param;
	double P; // / Pa
	double nu_div; // / (1 / m)
	double E_L, E_S; // / (W / m^2)
	i_0_z = n_eq + 1; // Index of the first value of altitude.
	i_0_param = i_0_z + n_eq + 1; // Index of the first parameter.
	nu_div = Y_0[i_0_param];
	for (int i = 0; i < n_eq; i++) {
		P = get_pressure(Y_0[i_0_z + i], Y_0[i]);
		E_L = longwave_irradiance(global_nu_min, nu_div, 55, n_eq, Y_0, Y_0 + i_0_z); // MC change to n_nu = 21 to have dnu ~ 100 / cm.
		E_S = shortwave_irradiance(nu_div, global_nu_max, 38, n_eq, Y_0, Y_0 + i_0_z); // MC change to n_nu = 978 to have dnu ~ 100 / cm.
		R[i] = - (E_L + E_S) / (Y_0[i_0_z + i + 1] - Y_0[i_0_z + i])
			/ (get_density(Y_0[i_0_z + i], Y_0[i]) * global_c_P_air);
	}
}
