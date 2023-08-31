#ifndef RADIATION_H
#define RADIATION_H

#include "constants.h"
#include "functions.h"
#include "configuration.h" // Include last to allow redefinitions.

struct absorber_t {
	int n_nu;
	double * nu; // / (1 / m)
	double * delta_nu; // / (1 / m)
};

double total_flux_longwave(double T, int n_nu, double nu[], double delta_nu[]);

double flux_longwave(double T, double nu, double delta_nu);

void delete_absorber(absorber_t absorber);

/*
Return wavenumber in 1 / m of intersection between spectral irradiances of Sun's surface at Earth's surface position and Earth's surface.
*/
double spectrum_division_nu();

/*
Return wavenumber in 1 / m of a photon in vacuum with energy Q / eV.
*/
double photon_wavenumber(double Q);

#endif /* RADIATION_H */
