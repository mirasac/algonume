#ifndef RADIATION_H
#define RADIATION_H

#include <cmath>

#include "constants.h"
#include "../../mclib/mclib.h"
#include "configuration.h" // Include last to allow redefinitions.

struct absorber_t {
	int n_nu;
	double * nu; // / (1 / m)
	double * delta_nu; // / (1 / m)
};

void delete_absorber(absorber_t absorber);

/*
Evaluate the spectral radiance of a blackbody for wavenumber nu / (1 / m) and temperature T / K.
Units: W m / (m^2 sr).
*/
double planck_law_nu(double nu, double T);

/*
Evaluate first derivative with respect to wavenumber nu / (1 / m) of the spectral radiance of a blackbody at temperature T / K.
Units: W m^2 / (m^2 sr).
*/
double planck_law_nu1(double nu, double T);

/*
Evaluate the spectral radiance of a blackbody for wavelenght lambda / m and temperature T / K.
Units: W / (m^2 m sr).
*/
double planck_law_lambda(double lambda, double T);

/*
Evaluate the average spectral radiance of a blackbody in the interval [nu, nu + dnu] and temperature T / K. Wavenumbers are expressed in unit 1 / m.
Units: W m / (m^2 sr).
*/
double planck_law_nu_average(double nu, double dnu, double T);

/*
Return wavenumber in 1 / m of intersection between spectral irradiances of Sun's surface at Earth's surface position and Earth's surface.
*/
double spectrum_division_nu();

/*
Return wavenumber in 1 / m of a photon in vacuum with energy Q / eV.
*/
double photon_wavenumber(double Q);

/*
Return longwave irradiance in W / m^2.
*/
double longwave_irradiance(double const nu_min, double const nu_max, int const n_nu, int const n_layers, double const T[], double const z[]);

#endif /* RADIATION_H */
