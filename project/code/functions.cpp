#include "functions.h"

double planck_law_nu(double nu, double T) {
	double factor = global_h * global_c / global_k_B;
	return 2.0 * global_h * global_c*global_c * nu*nu*nu / expm1(factor * nu / T);
}

double planck_law_lambda(double lambda, double T) {
	double factor = global_h * global_c / global_k_B;
	return 2.0 * global_h * global_c*global_c / lambda / lambda / lambda / (lambda*lambda * expm1(factor / (lambda * T)));
}

double planck_law_nu_average(double nu, double dnu, double T) {
	return gaussquad(planck_law_nu, nu, nu + dnu, GAUSSQUAD_INTERVALS, 2, T) / dnu;
}
