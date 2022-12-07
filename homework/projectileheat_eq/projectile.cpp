/*
MC things to do before PROD:
- Simplify fractions.
- Put meaningful comments and remove useless ones.
- Remove useless global constants.
- It seems that my bisection function has low accuracy, find why. A reason is the relatively low number of points x which I use for the integration. Another problem is the slow convergence to the correct results when incrementing the number of points x, probably due to a propagation of errors in bisection that may originate from the implementation of bisection which calls more times the functio, which is residual in this case, hence the number of internal calculations is high.
- Remove comments with units of measurement from code for PROD.
*/

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

// Do not copy into report.
#include "../../mclib/mclib.h"
//#define PROD
#ifdef N_PRECISION
#undef N_PRECISION
#endif /* N_PRECISION */

#define HOMEWORK_NAME "projectile"
#define N_PRECISION 12
#define TOLERANCE 1e-7

/* Global variables */
static double const global_g = 9.81;
static double const global_C = 0.1;
static double const global_v_0 = 18.0;
static double const global_x_0 = 0.0;
static double const global_L = 10.0;
static double const global_y_0 = 0.0;

/* Function prototypes */
void rhs(double const x, double const Y[], double R[]);
double residual(double phi);

int main() {
	// Setup.
	using namespace std;
	cout << setprecision(N_PRECISION) << scientific;
	int const n_eq = 3;
	int const n_x = 1000;  // MC number of points.
	double const dx = (global_L - global_x_0) / (n_x - 1);
	double x;
	double phi_0, phi_1;  // [rad]
	double Y_0[n_eq], Y_1[n_eq];
	ofstream plot_file;
	
	#ifndef PROD
	// Choose graphically the intervals of phi where phi_0 and phi_1 are searched.
	int const n_phi = 10;  // Number of intervals.
	double phi, dphi, phi_min, phi_max;  // [rad]
	double Y[n_eq];
	phi_min = TOLERANCE;
	//phi_max = (1.0) * M_PI / 2.0;  // Old maximum angle, there are NaNs.
	phi_max = (1.0 - TOLERANCE) * M_PI / 2.0;  // Old maximum angle.
	// MC from 65 there is divergence of tan and discontinuity of phi'=dphi/dx, it is when is reached phi < -pi/2. Find the value.
	// MC contnuare.
	//phi_max = 1023.0 / 1024.0 * M_PI / 2.0; //phi_max = 85 * M_PI / 180.0;
	dphi = (phi_max - phi_min) / n_phi;
	plot_file.open(HOMEWORK_NAME "_noprod_search.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "x y(x) v(x) phi(x)" << endl;
	for (int i_phi = 0; i_phi <= n_phi; i_phi++) {
		phi = phi_min + i_phi * dphi;
		// Initial values.
		Y[0] = global_y_0;  // y [m]
		Y[1] = global_v_0;  // v [m / s]
		Y[2] = phi;  // phi [rad]
		for (int i_x = 0 ; i_x < n_x; i_x++) {
			x = global_x_0 + i_x * dx;
			plot_file << x << ' ' << Y[0] << ' ' << Y[1] << ' ' << Y[2] << endl;
			rungekutta2(x, dx, Y, rhs, n_eq);
			//if (Y[2] <= -M_PI_2) break;  // Useful to avoid divergence to infinity of the trajectory.
		}
		plot_file << '\n' << endl;
	}
	plot_file.close();
	#endif /* PROD */
	
	// Find initial shooting angles.
	phi_0 = bisection(residual, 3.0 * M_PI / 30.0, 4.0 * M_PI / 30.0, TOLERANCE);
	phi_1 = bisection(residual, 8.0 * M_PI / 30.0, 9.0 * M_PI / 30.0, TOLERANCE);
	cout << "phi_0 = " << phi_0 / M_PI * 180.0 << " deg" << endl;
	cout << "phi_1 = " << phi_1 / M_PI * 180.0 << " deg" << endl;
	
	// Save points of trajectories to file.
	// Initial values.
	Y_0[0] = global_y_0;  // y [m]
	Y_0[1] = global_v_0;  // v [m / s]
	Y_0[2] = phi_0;  // phi [rad]
	Y_1[0] = global_y_0;  // y [m]
	Y_1[1] = global_v_0;  // v [m / s]
	Y_1[2] = phi_1;  // phi [rad]
	plot_file.open(HOMEWORK_NAME ".dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "x y_0(x) y_1(x)" << endl;
	for (int i_x = 0; i_x < n_x; i_x++) {
		x = global_x_0 + i_x * dx;
		plot_file << x << ' ' << Y_0[0] << ' ' << Y_1[0] << endl;
		rungekutta4(x, dx, Y_0, rhs, n_eq);
		rungekutta4(x, dx, Y_1, rhs, n_eq);
	}
	
	// Teardown.
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
	int const n_x = 100;
	double const dx = (global_L - global_x_0) / (n_x - 1);
	double x;
	double Y[n_eq];
	Y[0] = global_y_0;  // y [m]
	Y[1] = global_v_0;  // v [m / s]
	Y[2] = phi;  // phi [rad]
	for (int i_x = 0 ; i_x < n_x; i_x++) {
		x = global_x_0 + i_x * dx;
		rungekutta4(x, dx, Y, rhs, n_eq);
	}
	return Y[0];  // Boundary value is y(L) = 0.
}
