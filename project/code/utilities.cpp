#include "utilities.h"

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

double get_density(double z, double T) {
	return get_pressure(z, T) / global_R_m_air / T;
}
