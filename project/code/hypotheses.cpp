#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "constants.h"
#include "../../mclib/mclib.h"
#include "radiation.h"
#include "configuration.h" // Include last to allow redefinitions.

int main(int argc, char * argv[]) {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	
	// Check validity of EM spectrum.
	double M_sun, E_sun, M_earth, E_earth; // / (W / m^2)
	M_sun = global_sigma * global_T_sun*global_T_sun*global_T_sun*global_T_sun;
	E_sun = M_PI * gaussquad(spectral_irradiance_blackbody_nu, global_nu_min, global_nu_max, QUAD_INTERVALS, 2, global_T_sun);
	M_earth = global_sigma * global_T_earth*global_T_earth*global_T_earth*global_T_earth;
	E_earth = M_PI * gaussquad(spectral_irradiance_blackbody_nu, global_nu_min, global_nu_max, QUAD_INTERVALS, 2, global_T_earth);
	cout << "In EM spectrum [" << global_nu_min / 100.0 << " 1 / cm, " << global_nu_max / 100.0 << " 1 / cm]:" << endl;
	cout << "- ratio of Sun's surface calculated irradiance to theroretical radiant exitance is: " << E_sun / M_sun << endl;
	cout << "- ratio of Earth's surface calculated irradiance to theroretical radiant exitance is: " << E_earth / M_earth << endl;
	
	// Plot spectral irradiance of Sun's and Earth's surfaces as blackbodies.
	cout << endl;
	int const n_nu = 10000;
	double ratio;
	double nu, dnu; // / (1 / m)
	double E_sun_nu, E_earth_nu; // / (W m / m^2)
	ofstream file_plot;
	char filename_plot[] = DIR_DATA "/spectral_irradiance.dat";
	ratio = global_R_sun / global_au;
	dnu = (global_nu_max - global_nu_min) / n_nu;
	file_plot.open(filename_plot);
	file_plot << fixed << setprecision(N_PRECISION);
	file_plot << "#nu       E_sun(nu)    E_earth(nu)" << endl;
	file_plot << "#'1 / cm' 'W cm / m^2' 'W cm / m^2'" << endl;
	for (int i = 0; i <= n_nu; i++) {
		nu = global_nu_min + i * dnu;
		E_sun_nu = (1.0 - global_A) * ratio*ratio * M_PI * spectral_irradiance_blackbody_nu(nu, global_T_sun);
		E_earth_nu = M_PI * spectral_irradiance_blackbody_nu(nu, global_T_earth);
		nu /= 100.0; // Plot bandwidth in unit 1 / cm.
		E_sun_nu *= 100.0; // Plot irradiance in unit W cm / m^2.
		E_earth_nu *= 100.0; // Plot irradiance in unit W cm / m^2.
		file_plot << nu << ' ' << E_sun_nu << ' ' << E_earth_nu << '\n';
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
	E_sun_long = (1.0 - global_A) * ratio*ratio * M_PI * gaussquad(spectral_irradiance_blackbody_nu, global_nu_min, nu_div, QUAD_INTERVALS, 2, global_T_sun);
	E_earth_long = M_PI * gaussquad(spectral_irradiance_blackbody_nu, global_nu_min, nu_div, QUAD_INTERVALS, 2, global_T_earth);
	E_sun_short = (1.0 - global_A) * ratio*ratio * M_PI * gaussquad(spectral_irradiance_blackbody_nu, nu_div, global_nu_max, QUAD_INTERVALS, 2, global_T_sun);
	E_earth_short = M_PI * gaussquad(spectral_irradiance_blackbody_nu, nu_div, global_nu_max, QUAD_INTERVALS, 2, global_T_earth);
	cout << "In longwave bandwidth [" << global_nu_min / 100.0 << " 1 / cm, " << nu_div / 100.0 << " 1 / cm]:" << endl;
	cout << "- irradiances ratio of Sun's surface to total: " << E_sun_long / (E_sun_long + E_earth_long) << endl;
	cout << "- irradiances ratio of Sun's surface to whole spectrum Sun's surface: " << E_sun_long / (E_sun_long + E_sun_short) << endl;
	cout << "In shortwave bandwidth [" << nu_div / 100.0 << " 1 / cm, " << global_nu_max / 100.0 << " 1 / cm]:" << endl;
	cout << "- irradiances ratio of Earth's surface to total: " << E_earth_short / (E_sun_short + E_earth_short) << endl;
	cout << "- irradiances ratio of Earth's surface to whole spectrum Earth's surface: " << E_earth_short / (E_earth_long + E_earth_short) << endl;
	
	return 0;
}
