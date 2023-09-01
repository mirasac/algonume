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

#endif /* CONFIGURATION_H */
