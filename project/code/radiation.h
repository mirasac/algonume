#ifndef RADIATION_H
#define RADIATION_H

#include "constants.h"
#include "functions.h"

double total_flux_longwave(double T, int n_nu, double nu[], double dnu[]);

double flux_longwave(double T, double nu, double dnu);

#endif /* RADIATION_H */
