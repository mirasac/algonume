#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "constants.h"
#include "functions.h"
#include "../../mclib/mclib.h"

double spectral_irradiance_diff(double nu);

int main(int argc, char * argv[]) {
	using namespace std;
	cout << scientific << setprecision(N_PRECISION);
	
	// Plot spectral irradiance of Sun and Earth's surfaces as blackbodies.
	int const n_nu = 100000;
	double nu_sun, nu_earth, nu, dnu, nu_min, nu_max, nu_intersection; // [1 / m]
	double I_sun, I_earth; // [W / (m^2 m)]
	double ratio;
	ofstream plot_file;
	nu_sun = 2e6;
	nu_earth = 1e5;
	nu_min = 1e4;
	nu_max = 1e7;
	ratio = global_R_sun / global_au;
	dnu = (nu_max - nu_min) / n_nu;
	plot_file.open("spectral_irradiance.dat");
	plot_file << fixed << setprecision(N_PRECISION);
	plot_file << "#nu I_sun I_earth" << endl;
	for (int i = 0; i <= n_nu; i++) {
		nu = nu_min + i * dnu;
		I_sun = 0.1 * (1.0 - global_alpha) * ratio*ratio * M_PI * planck_law_nu(nu, global_T_sun);
		I_earth = M_PI * planck_law_nu(nu, global_T_earth);
		plot_file << nu / 1e2 << ' ' << ' ' << I_sun << ' ' << I_earth << '\n';
	}
	plot_file.close();
	
	nu_intersection = secant(spectral_irradiance_diff, 2e5, 3e5, 1e-12) / 100.0;
	cout << nu_intersection << endl; // MC debug.
	
	// Evaluate overlap of spectral irradiances.
	
	
	return 0;
}

double spectral_irradiance_diff(double nu) {
	double ratio = global_R_sun / global_au;
	return M_PI * (0.1 * (1.0 - global_alpha) * ratio*ratio * planck_law_nu(nu, global_T_sun) - planck_law_nu(nu, 288.0));
}
