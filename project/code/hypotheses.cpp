#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "configuration.h"
#include "constants.h"
#include "functions.h"
#include "../../mclib/mclib.h"

double spectral_irradiance_diff(double nu);
double spectral_irradiance_diff1(double nu);

int main(int argc, char * argv[]) {
	using namespace std;
	cout << scientific << setprecision(N_PRECISION);
	
	int const n_nu = 10000;
	double nu, dnu, nu_min, nu_max, nu_div; // / (1 / m)
	double I_sun, I_earth; // / (W m / m^2)
	double ratio;
	ofstream plot_file;
	nu_min = 1e4;
	nu_max = 1e7;
	ratio = global_R_sun / global_au;
	dnu = (nu_max - nu_min) / n_nu;
	
	// Plot spectral irradiance of Sun and Earth's surfaces as blackbodies.
	plot_file.open(DIR_DATA "/spectral_irradiance.dat");
	plot_file << fixed << setprecision(N_PRECISION);
	plot_file << "#nu I_sun I_earth" << endl;
	for (int i = 0; i <= n_nu; i++) {
		nu = nu_min + i * dnu;
		I_sun = (1.0 - global_alpha) * ratio*ratio * M_PI * planck_law_nu(nu, global_T_sun);
		I_earth = M_PI * planck_law_nu(nu, global_T_earth);
		nu /= 100.0; // Plot bandwidth in unit 1 / cm.
		I_sun *= 100.0; // Plot irradiance in unit W cm / m^2.
		I_earth *= 100.0; // Plot irradiance in unit W cm / m^2.
		plot_file << nu << ' ' << I_sun << ' ' << I_earth << '\n';
	}
	plot_file.close();
	
	// Find intersection between spectral irradiances.
	nu_div = newtonraphson(spectral_irradiance_diff, spectral_irradiance_diff1, 2e5, 3e5, 1e-10);
	cout << "Spectral irradiances intersect at nu_div = " << nu_div / 100.0 << " 1 / cm" << endl;
	
	// Evaluate overlap of spectral irradiances.
	
	
	return 0;
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
