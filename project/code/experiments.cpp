#include <cmath>
#include <iomanip>
#include <iostream>

#include "constants.h"
#include "../../mclib/mclib.h"
#include "radiation.h"
#include "configuration.h" // Include last to allow redefinitions.

int main(int argc, char * argv[]) {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	
	// Check validity of EM spectrum over a range of temperatures.
	int n_T = 1000;
	double T, dT, T_min, T_max; // / K
	double M, E; // / (W / m^2)
	T_min = 0.0;
	T_max = 1e3;
	dT = (T_max - T_min) / n_T;
	cout << scientific;
	cout << "#T   M(T)      E(T)      error" << endl;
	cout << "#'K' 'W / m^2' 'W / m^2' ''" << endl;
	for (int i = 0; i <= n_T; i++) {
		T = T_min + i * dT;
		M = global_sigma * T*T*T*T;
		E = M_PI * gaussquad(planck_law_nu, global_nu_min, global_nu_max, QUAD_INTERVALS, 2, T);
		cout << T << ' ' << M << ' ' << E << ' ' << fabs(E - M) / M << '\n';
	}
	cout << fixed;
	
	// Evaluate overlap of spectral irradiances for shortwave radiation.
	cout << endl;
	double ratio;
	double nu_div; // / (1 / m)
	double E_sun_short, E_sun_short_IR; // / (W / m^2)
	ratio = global_R_sun / global_au;
	nu_div = spectrum_division_nu();
	E_sun_short = (1.0 - global_A) * ratio*ratio * M_PI * gaussquad(planck_law_nu, nu_div, global_nu_max, QUAD_INTERVALS, 2, global_T_sun);
	E_sun_short_IR = (1.0 - global_A) * ratio*ratio * M_PI * gaussquad(planck_law_nu, nu_div, global_nu_IR_VIS, QUAD_INTERVALS, 2, global_T_sun);
	cout << "In shortwave bandwidth [" << nu_div / 100.0 << " 1 / cm, " << global_nu_max / 100.0 << " 1 / cm]:" << endl;
	cout << "- Sun's surface irradiances ratio of shortwave IR band to shortwave bandwidth: " << E_sun_short_IR / E_sun_short << endl;
	
	return 0;
}
