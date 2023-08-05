#include "functions.h"

double planck_law_nu(double nu, double T) {
	double factor = global_h * global_c * nu;
	return 2.0 * factor * global_c * nu*nu / expm1(factor / global_k_B / T);
}

double planck_law_nu_average(double nu, double dnu, double T) {
	return gaussquad2(planck_law_nu, nu, nu + dnu, GAUSS_QUAD_INTERVALS, GAUSS_QUAD_POINTS, T) / dnu;
}
