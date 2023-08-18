#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "configuration.h"
#include "constants.h"
#include "functions.h"
#include "../../mclib/mclib.h"

double spectral_irradiance_diff(double nu);

int main(int argc, char * argv[]) {
	using namespace std;
	cout << scientific << setprecision(N_PRECISION);
	
	// Plot spectral irradiance of Sun and Earth's surfaces as blackbodies.
	int const n_nu = 10000;
	double nu_sun, nu_earth, nu, dnu, nu_min, nu_max, nu_intersection; // / (1 / m)
	double I_sun, I_earth; // / (W m / m^2)
	double ratio;
	ofstream plot_file;
	nu_sun = 2e6;
	nu_earth = 1e5;
	nu_min = 1e4;
	nu_max = 1e7;
	ratio = global_R_sun / global_au;
	dnu = (nu_max - nu_min) / n_nu;
	plot_file.open(DIR_DATA "/spectral_irradiance.dat");
	plot_file << fixed << setprecision(N_PRECISION);
	plot_file << "#nu I_sun I_earth" << endl;
	for (int i = 0; i <= n_nu; i++) {
		nu = nu_min + i * dnu;
		I_sun = (1.0 - global_alpha) * ratio*ratio * M_PI * planck_law_nu(nu, global_T_sun);
		I_earth = M_PI * planck_law_nu(nu, global_T_earth);
		I_sun /= 32.0; // Scale for better plot comparison.
		nu /= 100.0; // Plot bandwidth in unit 1 / cm.
		plot_file << nu << ' ' << I_sun << ' ' << I_earth << '\n';
	}
	plot_file.close();
	
	nu_intersection = secant(spectral_irradiance_diff, 2e5, 3e5, 1e-12) / 100.0;
	cout << nu_intersection << endl; // MC debug.
	
	cout << (1.0 - global_alpha) * ratio*ratio * M_PI * gaussquad(planck_law_lambda, 1.0 / nu_max, 1.0 / nu_min, GAUSSQUAD_INTERVALS, GAUSSQUAD_POINTS, global_T_sun) << endl;
	cout << (1.0 - global_alpha) * ratio*ratio * M_PI * gaussquad(planck_law_nu, nu_min, nu_max, GAUSSQUAD_INTERVALS, GAUSSQUAD_POINTS, global_T_sun) << endl;
	
	// Evaluate overlap of spectral irradiances.
	
	
	return 0;
}

double spectral_irradiance_diff(double nu) {
	double ratio = global_R_sun / global_au;
	return M_PI * ((1.0 - global_alpha) * ratio*ratio * planck_law_nu(nu, global_T_sun) / 32.0 - planck_law_nu(nu, 288.0));
}
