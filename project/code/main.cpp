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
	double * z, * delta_z; // / m
	z = new double[global_N];
	delta_z = new double[global_N];
	set_layers_z_uniform(global_z_g, global_z_TOA, global_N, z, delta_z);
	
	// Configure absorbers.
	absorber_t H2O, CO2;
	H2O.nu = new double[]{14000.0, 160000.0, 376000.0, 535000.0, 725000.0};
	H2O.n_nu = 5;
	H2O.delta_nu = new double[H2O.n_nu]; // MC bandwidths missing.
	CO2.nu = new double[]{66700.0, 96000.0, 106000.0, 241000.0, 366000.0, 520000.0};
	CO2.n_nu = 6;
	CO2.delta_nu = new double[CO2.n_nu]; // MC bandwidths missing.
	
	// Grey atmosphere profile, plot analytical solution.
	int const n_tau = 10000;
	double F_0, D, tau, tau_min, tau_max, dtau, T4;
	ofstream plot_file;
	F_0 = (1.0 - global_alpha) * global_S * 0.25;
	D = 1.66;
	tau_min = 0.0;
	tau_max = 10.0;
	dtau = (tau_max - tau_min) / (n_tau - 1);
	plot_file.open("MC.dat");
	plot_file << fixed << setprecision(N_PRECISION);
	plot_file << "#tau T(tau)" << endl;
	for (int i = 0; i < n_tau; i++) {
		tau = tau_min + dtau * i;
		T4 = F_0 * 0.5 / global_sigma * (1.0 + D * tau);
		plot_file << tau << ' ' << pow(T4, 0.25) << '\n';
	}
	plot_file.close();
	
	int * n_absorbers;
	absorber_t ** absorbers;
	n_absorbers = new int[global_N];
	absorbers = new absorber_t*[global_N];
	set_absorbers_uniform(global_N, n_absorbers, absorbers, H2O, CO2);
	
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
	
	// Pre-evaluate reflectance, absorptance and internal transmittance.
	double * tau, * rho;
	tau = new double[global_N];
	rho = new double[global_N + 1];
	tau[0] = exp(-0.15); // MC test shortwave with constant optical depth delta = 0.15.
	rho[global_N] = global_A_g;
	for (int i = 0; i < global_N; i++) {
		tau[i] = tau[0];
		rho[i] = 0.0; // MC test shortwave with null reflectance.
	}
	
	// Initialise temperature output variable, set initial and boundary conditions.
	int n_T, i_0_z, i_0_param, i_0_tau, i_0_alpha;
	double T_TOA; // / K
	double * T; // / K
	i_0_z = global_N + 1; // Index of the first value of altitude.
	i_0_param = i_0_z + global_N + 1; // Index of the first parameter.
	i_0_tau = i_0_param + 1; // Index of the first value of internal transmittance.
	i_0_alpha = i_0_tau + global_N; // Index of the first value of absorptance shortwave.
	n_T = i_0_alpha + global_N; // Number of elements for the extended array of temperatures.
	T = new double[n_T]; // Use also for parameters.
	for (int i = 0; i < global_N; i++) {
		T[i] = global_T_earth;
		T[i_0_alpha + i] = absorptance_shortwave(i, global_N, tau, rho);
		T[i_0_z + i] = z[i];
	}
	T[global_N] = global_T_earth; // Temperature at ground level.
	T[i_0_z + global_N] = global_z_g;
	T[i_0_param] = nu_div;
	
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
	file_plot << "#'h' 'm' 'K' 'Pa' '1'   'K'";
	
	// Run model.
	i_t = 0;
	do {
		file_plot << '\n';
		t = i_t * dt;
		T_TOA = T[0];
		P_TOA = get_pressure(z[0], T_TOA);
		for (int i_z = 0; i_z < global_N; i_z++) {
			P = get_pressure(z[i_z], T[i_z]);
			file_plot << t << ' ' << z[i_z] << ' ' << T[i_z] << ' ' << P << ' ' << get_sigma(P, P_TOA) << ' ' << get_theta(T[i_z], P) << '\n';
		}
		file_plot << t << ' ' << global_z_g << ' ' << T[global_N] << ' ' << global_P_g << ' ' << get_sigma(global_P_g, P_TOA) << ' ' << get_theta(T[global_N], global_P_g) << '\n';
		eulerstep(t, dt, T, rhs, global_N);
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
	int i_0_z, i_0_param, i_0_tau, i_0_alpha;
	double const * z;
	double nu_div; // / (1 / m)
	double const * alpha;
	double E_L, E_S; // / (W / m^2)
	i_0_z = n_eq + 1; // Index of the first value of altitude.
	i_0_param = i_0_z + n_eq + 1; // Index of the first parameter.
	i_0_tau = i_0_param + 1; // Index of the first value of internal transmittance.
	i_0_alpha = i_0_tau + n_eq; // Index of the first value of absorptance shortwave.
	z = Y_0 + i_0_z;
	nu_div = Y_0[i_0_param];
	alpha = Y_0 + i_0_alpha;
	for (int i = 0; i < n_eq; i++) {
		E_L = irradiance_longwave(global_nu_min, nu_div, 55, n_eq, i, Y_0, z); // MC change to n_nu = 21 to have dnu ~ 100 / cm.
		E_S = global_S_0 * M_SQRT2 / 8.0 * alpha[i];
		R[i] = - (E_L + E_S) / (z[i + 1] - z[i]) / (get_density(z[i], Y_0[i]) * global_c_P_air);
	}
}
