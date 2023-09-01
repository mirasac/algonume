#ifndef UTILITIES_H
#define UTILITIES_H

#include <cmath>

#include "constants.h"
#include "radiation.h"

void set_layers_z_uniform(double z_min, double z_max, int n_layers, double z[], double delta_z[]);

void set_absorbers_uniform(int n_layers, int n_absorbers[], absorber_t * absorbers[], absorber_t a1, absorber_t a2);

double get_pressure(double z, double T);

double get_altitude(double P, double T);

double get_sigma(double P, double P_TOA);

double get_theta(double T, double P);

double get_density(double z, double T);

#endif /* UTILITIES_H */
