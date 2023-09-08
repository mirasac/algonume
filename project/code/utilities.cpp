#include "utilities.h"

void set_layers_uniform(double x_min, double x_max, int n_layer, double x[], double dx[]) {
	dx[0] = (x_max - x_min) / n_layer;
	x[0] = x_max;
	for (int i = 1; i < n_layer; i++) {
		dx[i] = dx[0];
		x[i] = x[0] - dx[0] * i;
	}
}

void set_layers_z_uniform(double z_min, double z_max, int n_layer, layer_t layer[]) {
	layer[0].dz = (z_max - z_min) / n_layer;
	layer[0].z = z_max;
	for (int i = 1; i < n_layer; i++) {
		layer[i].dz = layer[0].dz;
		layer[i].z = layer[0].z - layer[0].dz * i;
	}
}

void set_absorbers_uniform(int n_layer, int n_absorbers[], absorber_t * absorbers[], absorber_t a1, absorber_t a2) {
	n_absorbers[0] = 2;
	absorbers[0] = new absorber_t[n_layer * n_absorbers[0]];
	absorbers[0][0] = a1;
	absorbers[0][1] = a2;
	for (int i = 1; i < n_layer; i++) {
		n_absorbers[i] = n_absorbers[0];
		absorbers[i] = absorbers[0] + i * n_absorbers[0];
		absorbers[i][0] = a1;
		absorbers[i][1] = a2;
	}
}

double get_pressure(double z) {
	return global_P_g * exp(- (z - global_z_g) / global_z_0);
}

double get_altitude(double P) {
	return global_z_g - global_z_0 * log(P / global_P_g);
}

double get_altitude(double P, double T) {
	return global_z_g - global_R_m * T / global_g * log(P / global_P_g);
}

double get_sigma(double P, double P_TOA) {
	return (P - P_TOA) / (global_P_g - P_TOA);
}

double get_theta(double T, double P) {
	return T * pow(global_P_0 / P, global_R_m / global_c_P_air);
}

double get_density(double z, double T) {
	return get_pressure(z) / global_R_m / T;
}

double get_optical_depth_P(double P, double P_TOA) {
	return global_mu_m / global_g * (P - P_TOA);
}

double get_optical_depth_z(double z, double P_TOA) {
	return get_optical_depth_P(get_pressure(z), P_TOA);
}
