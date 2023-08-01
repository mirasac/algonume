#include "radiation.h"

double total_flux_longwave(int n_nu, double nu[], double dnu[], double T){
	double _return = 0.0;
	for (int i = 0; i < n_nu; i++) {
		_return += flux_longwave(T, nu[i], dnu[i]) * dnu[i];
	}
	return _return;
}

double flux_longwave(double nu, double dnu, double T) {
	double _return; // MC continue with parametrization of T(z, z').
	return _return;
}
