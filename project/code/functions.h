#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>

#include "constants.h"
#include "../../mclib/mclib.h"

#define GAUSSQUAD_POINTS 2
#define GAUSSQUAD_INTERVALS 100000

/*
Evaluate the spectral radiance of a blackbody for wavenumber nu [1 / m] and temperature T [K].
Units: W m / (m^2 sr).
*/
double planck_law_nu(double nu, double T);

/*
Evaluate the spectral radiance of a blackbody for wavelenght lambda [m] and temperature T [K].
Units: W / (m^2 m sr).
*/
double planck_law_lambda(double lambda, double T);

/*
Evaluate the average spectral radiance of a blackbody in the interval [nu, nu + dnu] and temperature T [K]. Wavenumbers are expressed in unit [1 / m].
Units: W m / (m^2 sr).
*/
double planck_law_nu_average(double nu, double dnu, double T);

#endif /* FUNCTIONS_H */
