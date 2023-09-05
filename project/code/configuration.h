#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#define DIR_DATA "../data"

#ifdef N_PRECISION
#undef N_PRECISION
#endif /* N_PRECISION */
#define N_PRECISION 6

#define TOLERANCE 1e-6

#define GAUSSQUAD_POINTS 2
#define QUAD_INTERVALS 100000

// Arbitrary values.
double const global_delta_g = 0.8;
int const global_N = 100;
double const global_nu_min = 1e4; // 1 / m
double const global_nu_max = 1e7; // 1 / m
double const global_P_0 = 1e5; // Pa
double const global_P_TOA = 3.0; // Pa
double const global_z_0 = 2000.0; // / m
double const global_z_TOA = 55000.0; // m

double const global_mu_m = global_delta_g * global_g / (global_P_g - global_P_TOA);

#endif /* CONFIGURATION_H */
