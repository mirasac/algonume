#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "constants.h"

#define DIR_DATA "../data"

#ifdef N_PRECISION
#undef N_PRECISION
#endif /* N_PRECISION */
#define N_PRECISION 6

#define TOLERANCE 1e-6

#define GAUSSQUAD_POINTS 2
#define QUAD_INTERVALS 100000

#define N_STABILITY 25

// Arbitrary values.
int const global_N = 100;
double const global_nu_min = 1e4; // 1 / m
double const global_nu_max = 1e7; // 1 / m
double const global_P_0 = 1e5; // Pa
double const global_P_TOA = 3.0; // Pa
double const global_S_t = 0.25 * (1.0 - const_A) * const_S_0; // / (W / m^2)
double const global_delta_g = 2.0 * (const_sigma / global_S_t * const_T_g*const_T_g*const_T_g*const_T_g - 1.0) / const_D;
double const global_mu_m = global_delta_g * const_g / (const_P_g - global_P_TOA); // / (m^2 / kg)
double const global_z_0 = 2000.0; // m

#endif /* CONFIGURATION_H */
