/*
MC things to do before PROD:
- Simplify fractions.
- Put meaningful comments and remove useless ones.
- Remove useless global constants.
- It seems that my bisection function has low accuracy, find why. A reason is the relatively low number of points x which I use for the integration. Another problem is the slow convergence to the correct results when incrementing the number of points x, probably due to a propagation of errors in bisection that may originate from the implementation of bisection which calls more times the functio, which is residual in this case, hence the number of internal calculations is high.
*/

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../../mclib/mclib.h"

#ifdef N_PRECISION
#undef N_PRECISION
#define N_PRECISION 7
#endif /* N_PRECISION */

#define PROD

static double const global_g = 9.81;  // [m / s^2]
static double const global_C = 0.1;  // [m^-1]
static double const global_v_0 = 18.0;  // [m / s]
static double const global_x_min = 0.0;  // [m]
static double const global_L = 10.0;  // [m]
static double const global_y_0 = 0.0;  // [m]
//static double const global_y_L = 0.0;  // [m] MC not needed.

void rhs(double const x, double const Y[], double R[]);

double residual(double phi);

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	int const n_eq = 3;
	int n_x;
	double const tolerance = 1e-7;
	double phi_0, phi_1;  // [rad]
	double x_min, x_max, dx, x;  // [m]
	double Y_0[n_eq], Y_1[n_eq];
	ofstream plot_file;
	n_x = 100;  // Number of points.
	x_min = global_x_min;
	x_max = global_L;
	dx = (x_max - x_min) / (n_x - 1);
	#ifndef PROD
	// Choose graphically the intervals of phi where phi_0 and phi_1 are searched.
	int n_phi;
	double phi, dphi, phi_max, phi_min, phi_prev;  // [rad]
	double Y[n_eq];
	n_phi = 10;  // Number of intervals.
	phi_min = 0.0;
	phi_max = M_PI / 3.0;  // MC from 65 there is divergence of tan and discontinuity of phi'=dphi/dx, pi / 2 not chosen beacuse it gives inf.
	dphi = (phi_max - phi_min) / n_phi;
	plot_file.open("projectile_noprod_1.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "x y(x)" << endl;
	for (int i_phi = 0; i_phi <= n_phi; i_phi++) {
		phi = phi_min + i_phi * dphi;
		// Initial values.
		Y[0] = global_y_0;  // y [m]
		Y[1] = global_v_0;  // v [m / s]
		Y[2] = phi;  // phi [rad]
		for (int i_x = 0 ; i_x < n_x; i_x++) {
			x = x_min + i_x * dx;
			plot_file << x << ' ' << Y[0] << endl;
			phi_prev = Y[2];
			rungekutta4(x, dx, Y, rhs, n_eq);
			//if (phi_prev * Y[2] < 0.0) break;  // Useful to avoid divergence to infinity of the trajectory.
		}
		plot_file << '\n' << endl;
	}
	plot_file.close();
	#endif /* PROD */
	phi_0 = bisection(residual, 3.0 * M_PI / 30.0, 4.0 * M_PI / 30.0, tolerance);
	phi_1 = bisection(residual, 8.0 * M_PI / 30.0, 9.0 * M_PI / 30.0, tolerance);
	cout << "phi_0 [deg] = " << phi_0 / M_PI * 180.0 << endl;
	cout << "phi_1 [deg] = " << phi_1 / M_PI * 180.0 << endl;
	Y_0[0] = global_y_0;  // y [m]
	Y_0[1] = global_v_0;  // v [m / s]
	Y_0[2] = phi_0;  // phi [rad]
	Y_1[0] = global_y_0;  // y [m]
	Y_1[1] = global_v_0;  // v [m / s]
	Y_1[2] = phi_1;  // phi [rad]
	plot_file.open("projectile.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "x y_0(x) y_1(x)" << endl;
	for (int i_x = 0; i_x < n_x; i_x++) {
		x = x_min + i_x * dx;
		plot_file << x << ' ' << Y_0[0] << ' ' << Y_1[0] << endl;
		rungekutta4(x, dx, Y_0, rhs, n_eq);
		rungekutta4(x, dx, Y_1, rhs, n_eq);
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
	int const n_eq = 3;
	int n_x;
	double x_min, x_max, x, dx;  // [m]
	double Y[n_eq];
	n_x = 100;
	x_min = global_x_min;
	x_max = global_L;
	dx = (x_max - x_min) / (n_x - 1);
	Y[0] = global_y_0;  // y [m]
	Y[1] = global_v_0;  // v [m / s]
	Y[2] = phi;  // phi [rad]
	for (int i_x = 0 ; i_x < n_x; i_x++) {
		x = x_min + i_x * dx;
		rungekutta4(x, dx, Y, rhs, n_eq);
	}
	return Y[0];  // Boundary condition is y(L) = 0.
}
