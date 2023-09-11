#ifndef UTILITIES_H
#define UTILITIES_H

#include <cmath>

#include "constants.h"
#include "radiation.h"
#include "configuration.h" // Include last to allow redefinitions.

struct layer_t {
	double z; // / m
	double dz; // / m
	double T; // / K
	int n_gas;
	double * absorber_t;
};

void set_layers_uniform(double x_min, double x_max, int n_layer, double x[], double dx[]);

void set_layers_z_uniform(double z_min, double z_max, int n_layer, layer_t layer[]);

void set_absorbers_uniform(int n_layer, int n_absorbers[], absorber_t * absorbers[], absorber_t a1, absorber_t a2);

double get_pressure(double z);

/*
Pressure P / Pa must be positive.
*/
double get_altitude(double P);

double get_altitude(double P, double T);

double get_sigma(double P, double P_TOA);

/*
Potential temperature, proportional to temperature T / K and a negative power of pressure P / Pa.
*/
double get_theta(double T, double P);

double get_density(double z, double T);

double get_optical_depth_P(double P, double P_TOA);

double get_optical_depth_z(double z, double P_TOA);

/*
Return analytical solution for normalised temperature in radiative equilibrium.
*/
double temperature_norm(double delta);

/*
Return analytical solution for normalised upward irradiance in radiative equilibrium.
*/
double irradiance_upward_norm(double delta);

/*
Return analytical solution for normalised downward irradiance in radiative equilibrium.
*/
double irradiance_downward_norm(double delta);

#endif /* UTILITIES_H */
