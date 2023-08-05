#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>

#include "constants.h"
#include "../../mclib/mclib.h"

#define GAUSS_QUAD_POINTS 4
#define GAUSS_QUAD_INTERVALS 10

/*
Evaluate the spectral radiance of a blackbody for wavenumber nu [1 / m] and temperature T [K].
*/
double planck_law_nu(double nu, double T);

/*
Evaluate the average spectral radiance of a blackbody in the interval [nu, nu + dnu] and temperature T [K]. Wavenumbers are expressed in unit [1 / m].
*/
double planck_law_nu_average(double nu, double dnu, double T);

#endif /* FUNCTIONS_H */
