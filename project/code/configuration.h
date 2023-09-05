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
double const global_D = 1.66;
int const global_N = 20;
double const global_nu_min = 1e4; // 1 / m
double const global_nu_max = 1e7; // 1 / m
double const global_P_TOA = 100.0; // Pa
double const global_z_TOA = 55000.0; // m

#endif /* CONFIGURATION_H */
