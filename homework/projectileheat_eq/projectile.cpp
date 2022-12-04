#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../../mclib/mclib.h"

#ifdef N_PRECISION
#undef N_PRECISION
#define N_PRECISION 7
#endif /* N_PRECISION */

static double const global_g = 9.81;  // [m / s^2]
static double const global_C = 0.1;  // [m^-1]

void rhs(double const x, double const Y[], double R[]);

double residual(double phi);

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	int const n_eq = 3;
	int n_phi, n_x;
	double const v_0 = 18;  // [m / s]
	double const L = 10;  // [m]
	double phi, dphi, phi_min, phi_max;  // [rad]
	double x_min, x_max, dx, x;  // [m]
	double Y[n_eq];
	ofstream plot_file;
	n_phi = 90;
	phi_max = 99.0 / 100.0 * (M_PI / 2.0);
	phi_min = 1.0 / 100.0;
	dphi = (phi_max - phi_min) / n_phi;
	n_x = 100;
	dx = (x_max - x_min) / n_x;
	for (int i_phi = 0; i_phi < n_phi; i_phi++) {
		phi = phi_min + i_phi * dphi;
		// MC start copy into residual.
		// Initial values.
		Y[0] = 0.0;  // y [m]
		Y[1] = v_0;  // v [m / s]
		Y[2] = phi;  // phi [degrees]
		plot_file.open("projectile.dat");
		plot_file << setprecision(N_PRECISION) << scientific;
		plot_file << "x y(x)" << endl; 
		for (int i_x = 0 ; i_x < n_x; i_x++) {
			x = x_min + i_x * dx;
			plot_file << x << ' ' << Y[0] << endl;
			rungekutta4(x, dx, Y, rhs, n_eq);
		}
		// MC end copy into residual.
		plot_file << '\n' << endl;
	}
	plot_file.close();
	return 0;
}

void rhs(double const x, double const Y[], double R[]) {
	R[0] = tan(Y[2]);
	R[1] = - global_g / Y[1] * tan(Y[2]) - global_C * Y[1] / cos(Y[2]);
	R[2] = - global_g / (Y[1]*Y[1]);
}

double residual(double phi) {
	// MC copy here.
	return 0.0;
}
