#ifndef RADIATION_H
#define RADIATION_H

#include "constants.h"
#include "functions.h"

struct absorber_t {
	int n_nu;
	double * nu; // / (1 / cm)
	double * delta_nu; // / (1 / cm)
};

double total_flux_longwave(double T, int n_nu, double nu[], double delta_nu[]);

double flux_longwave(double T, double nu, double delta_nu);

void delete_absorber(absorber_t absorber);

#endif /* RADIATION_H */
