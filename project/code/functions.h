#ifndef UTILITIES_H
#define UTILITIES_H

#include <cmath>

#include "constants.h"
#include "../../mclib/mclib.h"

#define GAUSS_QUAD_POINTS 4
#define GAUSS_QUAD_INTERVALS 10

double planck_law(double nu, double T);

double planck_law_average(double nu, double dnu, double T);

#endif /* UTILITIES_H */
