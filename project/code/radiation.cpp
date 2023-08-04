#include "radiation.h"

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
