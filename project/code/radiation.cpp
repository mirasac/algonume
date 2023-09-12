#include "radiation.h"

void delete_absorber(absorber_t absorber) {
	delete[] absorber.nu;
	delete[] absorber.dnu;
}

double spectral_irradiance_blackbody_nu(double nu, double T) {
	double factor = const_h * const_c / const_k_B;
	return 2.0 * const_h * const_c*const_c * nu*nu*nu / expm1(factor * nu / T);
}

double spectral_irradiance_blackbody_nu1(double nu, double T) {
	double factor = const_h * const_c / const_k_B;
	return spectral_irradiance_blackbody_nu(nu, T) / nu * (3.0 + factor * nu / T / expm1(-factor * nu / T));
}

double spectral_irradiance_blackbody_lambda(double lambda, double T) {
	double factor = const_h * const_c / const_k_B;
	return 2.0 * const_h * const_c*const_c / lambda / lambda / lambda / (lambda*lambda * expm1(factor / (lambda * T)));
}

double irradiance_blackbody_average_nu(double nu, double dnu, double T) {
	return gaussquad(spectral_irradiance_blackbody_nu, nu, nu + dnu, QUAD_INTERVALS, 2, T) / dnu;
}

double spectral_irradiance_diff(double nu) {
	double ratio = const_R_sun / const_au;
	return M_PI * ((1.0 - const_A) * ratio*ratio * spectral_irradiance_blackbody_nu(nu, const_T_sun) - spectral_irradiance_blackbody_nu(nu, const_T_g));
}

double spectral_irradiance_diff1(double nu) {
	double ratio, factor;
	ratio = const_R_sun / const_au;
	factor = const_h * const_c / const_k_B;
	return M_PI * (3.0 / nu * spectral_irradiance_diff(nu) + factor * (
		(1.0 - const_A) * ratio*ratio
		* spectral_irradiance_blackbody_nu(nu, const_T_sun) / (const_T_sun * expm1(-factor * nu / const_T_sun))
		- spectral_irradiance_blackbody_nu(nu, const_T_g) / (const_T_g * expm1(-factor * nu / const_T_g))
	));
}

double spectrum_division_nu() {
	return newtonraphson(spectral_irradiance_diff, spectral_irradiance_diff1, 2e5, 3e5, TOLERANCE);
}

double wavenumber_photon(double Q) {
	return const_h * const_c / (const_e * Q);
}

// MC test with constant optical depth delta = 0.15.
double spectral_longwave_irradiance(double nu, double T_g) {
	return M_PI * spectral_irradiance_blackbody_nu(nu, T_g) * exp(-0.15);
}

// MC beware that here I am supposing that all layers have same thickness.
double irradiance_longwave(double const nu_min, double const nu_max, int const n_nu, int const i, int const n_layer, double const T[], double const z[]) {
	return trapezioidalquad(spectral_longwave_irradiance, nu_min, nu_max, n_nu, T[n_layer]);
}

double absorptance_shortwave(int const i, int const n_layer, double const tau[], double const rho[]) {
	double p1, p2, s1;
	p1 = 1.0;
	for (int j = 0; j <= i; j++) {
		p1 *= tau[j] * (1.0 - rho[j]);
	}
	s1 = 0.0;
	for (int n = i; n <= n_layer - 1; n++) {
		p2 = rho[n + 1];
		for (int j = i + 1; j <= n; j++) {
			p2 *= tau[j]*tau[j] * (1.0 - rho[j]);
		}
	}
	return (1.0 - tau[i]) * p1 * (1.0 / tau[i] + s1);
}

//nu_div, global_nu_max, 38, n_eq, Y_0, Y_0 + i_0_z); // MC change to n_nu = 978 to have dnu ~ 100 / cm.
