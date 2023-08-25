#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "constants.h"
#include "functions.h"
#include "../../mclib/mclib.h"
#include "radiation.h"
#include "configuration.h" // Load last to redefine some things.

int main(int argc, char * argv[]) {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	
	// Plot spectral irradiance of Sun and Earth's surfaces as blackbodies.
	int const n_nu = 10000;
	double nu, dnu, nu_min, nu_max; // / (1 / m)
	double E_sun, E_earth; // / (W m / m^2)
	double ratio;
	ofstream file_plot;
	char filename_plot[] = DIR_DATA "/spectral_irradiance.dat";
	nu_min = 1e4;
	nu_max = 1e7;
	ratio = global_R_sun / global_au;
	dnu = (nu_max - nu_min) / n_nu;
	file_plot.open(filename_plot);
	file_plot << fixed << setprecision(N_PRECISION);
	file_plot << "#nu E_sun E_earth" << endl;
	for (int i = 0; i <= n_nu; i++) {
		nu = nu_min + i * dnu;
		E_sun = (1.0 - global_alpha) * ratio*ratio * M_PI * planck_law_nu(nu, global_T_sun);
		E_earth = M_PI * planck_law_nu(nu, global_T_earth);
		nu /= 100.0; // Plot bandwidth in unit 1 / cm.
		E_sun *= 100.0; // Plot irradiance in unit W cm / m^2.
		E_earth *= 100.0; // Plot irradiance in unit W cm / m^2.
		file_plot << nu << ' ' << E_sun << ' ' << E_earth << '\n';
	}
	file_plot.close();
	cout << "Spectral irradiances plotted in file " << filename_plot << endl;
	
	// Find intersection between spectral irradiances.
	cout << endl;
	double nu_div; // / (1 / m)
	nu_div = spectrum_division_nu();
	cout << "Spectral irradiances intersect at nu_div = " << nu_div / 100.0 << " 1 / cm" << endl;
	
	// Evaluate overlap of spectral irradiances.
	cout << endl;
	double E_sun_long, E_earth_long, E_sun_short, E_earth_short; // / (W / m^2)
	E_sun_long = ratio*ratio * M_PI * gaussquad(planck_law_nu, nu_min, nu_div, GAUSSQUAD_INTERVALS, 2, global_T_sun);
	E_earth_long = M_PI * gaussquad(planck_law_nu, nu_min, nu_div, GAUSSQUAD_INTERVALS, 2, global_T_earth);
	E_sun_short = ratio*ratio * M_PI * gaussquad(planck_law_nu, nu_div, nu_max, GAUSSQUAD_INTERVALS, 2, global_T_sun);
	E_earth_short = M_PI * gaussquad(planck_law_nu, nu_div, nu_max, GAUSSQUAD_INTERVALS, 2, global_T_earth);
	cout << "In bandwidth [" << nu_min / 100.0 << " 1 / cm, " << nu_div / 100.0 << " 1 / cm]:" << endl;
	cout << " - irradiances ratio of Sun's surface to total: " << E_sun_long / (E_sun_long + E_earth_long) << endl;
	cout << " - irradiances ratio of Sun's surface to whole spectrum Sun's surface: " << E_sun_long / (E_sun_long + E_sun_short) << endl;
	cout << "In bandwidth [" << nu_div / 100.0 << " 1 / cm, " << nu_max / 100.0 << " 1 / cm]:" << endl;
	cout << " - irradiances ratio of Earth's surface to total: " << E_earth_short / (E_sun_short + E_earth_short) << endl;
	cout << " - irradiances ratio of Earth's surface to whole spectrum Earth's surface: " << E_earth_short / (E_earth_long + E_earth_short) << endl;
	
	return 0;
}
