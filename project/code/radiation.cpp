#include "radiation.h"

double planck_law_nu(double nu, double T) {
	double factor = global_h * global_c / global_k_B;
	return 2.0 * global_h * global_c*global_c * nu*nu*nu / expm1(factor * nu / T);
}

double planck_law_nu1(double nu, double T) {
	double factor = global_h * global_c / global_k_B;
	return planck_law_nu(nu, T) / nu * (3.0 + factor * nu / T / expm1(-factor * nu / T));
}

double planck_law_lambda(double lambda, double T) {
	double factor = global_h * global_c / global_k_B;
	return 2.0 * global_h * global_c*global_c / lambda / lambda / lambda / (lambda*lambda * expm1(factor / (lambda * T)));
}

double planck_law_nu_average(double nu, double dnu, double T) {
	return gaussquad(planck_law_nu, nu, nu + dnu, QUAD_INTERVALS, 2, T) / dnu;
}

double total_flux_longwave(int n_nu, double nu[], double delta_nu[], double T){
	double _return = 0.0;
	for (int i = 0; i < n_nu; i++) {
		_return += flux_longwave(T, nu[i], delta_nu[i]) * delta_nu[i];
	}
	return _return;
}

double flux_longwave(double nu, double delta_nu, double T) {
	double _return; // MC continue with parametrization of T(z, z').
	return _return;
}

void delete_absorber(absorber_t absorber) {
	delete[] absorber.nu;
	delete[] absorber.delta_nu;
}

double spectral_irradiance_diff(double nu) {
	double ratio = global_R_sun / global_au;
	return M_PI * ((1.0 - global_alpha) * ratio*ratio * planck_law_nu(nu, global_T_sun) - planck_law_nu(nu, global_T_earth));
}

double spectral_irradiance_diff1(double nu) {
	double ratio, factor;
	ratio = global_R_sun / global_au;
	factor = global_h * global_c / global_k_B;
	return M_PI * (3.0 / nu * spectral_irradiance_diff(nu) + factor * (
		(1.0 - global_alpha) * ratio*ratio
		* planck_law_nu(nu, global_T_sun) / (global_T_sun * expm1(-factor * nu / global_T_sun))
		- planck_law_nu(nu, global_T_earth) / (global_T_earth * expm1(-factor * nu / global_T_earth))
	));
}

double spectrum_division_nu() {
	return newtonraphson(spectral_irradiance_diff, spectral_irradiance_diff1, 2e5, 3e5, TOLERANCE);
}

double photon_wavenumber(double Q) {
	return global_h * global_c / (global_e * Q);
}
