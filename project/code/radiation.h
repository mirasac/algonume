#ifndef RADIATION_H
#define RADIATION_H

#include <cmath>

#include "constants.h"
#include "../../mclib/mclib.h"
#include "configuration.h" // Include last to allow redefinitions.

struct absorber_t {
	int n_nu;
	double * nu; // / (1 / m)
	double * dnu; // / (1 / m)
};

void delete_absorber(absorber_t absorber);

/*
Evaluate the spectral radiance of a blackbody for wavenumber nu / (1 / m) and temperature T / K.
Units: W m / (m^2 sr).
*/
double spectral_irradiance_blackbody_nu(double nu, double T);

/*
Evaluate first derivative with respect to wavenumber nu / (1 / m) of the spectral radiance of a blackbody at temperature T / K.
Units: W m^2 / (m^2 sr).
*/
double spectral_irradiance_blackbody_nu1(double nu, double T);

/*
Evaluate the spectral radiance of a blackbody for wavelenght lambda / m and temperature T / K.
Units: W / (m^2 m sr).
*/
double spectral_irradiance_blackbody_lambda(double lambda, double T);

/*
Evaluate the average spectral radiance of a blackbody in the interval [nu, nu + dnu] and temperature T / K. Wavenumbers are expressed in unit 1 / m.
Units: W m / (m^2 sr).
*/
double irradiance_blackbody_average_nu(double nu, double dnu, double T);

/*
Return wavenumber in 1 / m of intersection between spectral irradiances of Sun's surface at Earth's surface position and Earth's surface.
*/
double spectrum_division_nu();

/*
Return wavenumber in 1 / m of a photon in vacuum with energy Q / eV.
*/
double wavenumber_photon(double Q);

/*
Return longwave irradiance in W / m^2.
*/
double irradiance_longwave(double const nu_min, double const nu_max, int const n_nu, int const i, int const n_layer, double const T[], double const z[]);

/*
Return total absorptance due to shortwave radiation.

Array rho have length n_layer + 1, the element in position n_layer is the surface albedo.
*/
double absorptance_shortwave(int const i, int const n_layer, double const tau[], double const rho[]);

#endif /* RADIATION_H */
