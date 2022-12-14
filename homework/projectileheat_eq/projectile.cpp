#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

// Do not copy into report.
#include "../../mclib/mclib.h"
#define PROD
#ifdef N_PRECISION
#undef N_PRECISION
#endif /* N_PRECISION */

#define HOMEWORK_NAME "projectile"
#define TOLERANCE 1e-7
#define N_PRECISION 9

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
	int const n_x = 1025;
	double const dx = (global_L - global_x_0) / (n_x - 1);
	double x;
	double phi_0, phi_1;  // Unit of measurement: radian.
	double Y_0[n_eq], Y_1[n_eq];
	ofstream plot_file;
	
	#ifndef PROD
	// Choose graphically the intervals of phi where phi_0 and phi_1 are searched.
	int n_phi = 8;  // Number of intervals.
	double phi, dphi, phi_min, phi_max;  // [rad]
	double Y[n_eq];
	phi_min = TOLERANCE;
	phi_max = (1.0 - TOLERANCE) * M_PI / 2.0;  // Old maximum angle, divergence of trajectory arises due to numerical procedure.
	// MC from 65 there is divergence of tan and discontinuity of phi'=dphi/dx, it is when is reached phi < -pi/2. Find the value.
	phi_max = M_PI / 3.0;  // New maximum angle, before divergence is displayed.
	dphi = (phi_max - phi_min) / n_phi;
	plot_file.open(HOMEWORK_NAME "_noprod_search.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	//plot_file << "x y(x) v(x) phi(x)" << endl;  // Commented out because not good for gnuplot's index command.
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
			//if (Y[2] <= -M_PI_2) break;  // Useful to avoid plotting the divergence.
		}
		plot_file << '\n' << endl;
	}
	plot_file.close();

	// Show plot of residuals.
	n_phi = 1024;
	dphi = (phi_max - phi_min) / n_phi;
	plot_file.open(HOMEWORK_NAME "_noprod_residual.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "x y(x) v(x) phi(x)" << endl;
	for (int i_phi = 0; i_phi <= n_phi; i_phi++) {
		phi = phi_min + i_phi * dphi;
		plot_file << phi << ' ' << residual(phi) << endl;
	}
	plot_file.close();
	#endif /* PROD */
	
	// Find initial shooting angles.
	phi_0 = bisection(residual, M_PI / 12.0, M_PI / 8.0, TOLERANCE);
	phi_1 = bisection(residual, 7.0 * M_PI / 24.0, M_PI / 3.0, TOLERANCE);
	cout << "phi_0 = " << 180.0 / M_PI * phi_0 << " deg" << endl;
	cout << "phi_1 = " << 180.0 / M_PI * phi_1 << " deg" << endl;

	#ifndef PROD
	// Check relative error.
	cout << fabs(180.0 / M_PI * (phi_0 + phi_1) - 90) / TOLERANCE << endl;

	// Check using bracketing.
	int const n = 128;
	int n_bracketing;
	double phi_L[n], phi_R[n];
	n_bracketing = bracketing(residual, phi_min, phi_max, n, phi_L, phi_R);
	for (int i = 0; i < n_bracketing; i++) {
		cout << '[' << phi_L[i] << " rad, " << phi_R[i] << " rad]: " << 180.0 / M_PI * bisection(residual, phi_L[i], phi_R[i], TOLERANCE) << " deg" << endl;
	}
	#endif /* PROD */
	
	// Evaluate trajectories and save points to file.
	Y_0[0] = global_y_0;  // Physical quantity y [m]
	Y_0[1] = global_v_0;  // Physical quantity v [m / s]
	Y_0[2] = phi_0;  // Physical quantity phi [rad]
	Y_1[0] = global_y_0;  // Physical quantity y [m]
	Y_1[1] = global_v_0;  // Physical quantity v [m / s]
	Y_1[2] = phi_1;  // Physical quantity phi [rad]
	plot_file.open(HOMEWORK_NAME ".dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "x y_0(x) y_1(x)" << endl;
	for (int i_x = 0; i_x < n_x; i_x++) {
		x = global_x_0 + i_x * dx;
		plot_file << x << ' ' << Y_0[0] << ' ' << Y_1[0] << endl;
		rungekutta4(x, dx, Y_0, rhs, n_eq);
		rungekutta4(x, dx, Y_1, rhs, n_eq);
	}

	#ifndef PROD
	// Analytical results for initial shooting angles in free fall.
	double const tmp_arg = global_g * global_L / (global_v_0*global_v_0);
	double phi_rad;
	cout << "Analytical results" << endl;
	phi_rad = 0.5 * asin(tmp_arg);
	cout << "arcsin phi_0 " << phi_rad << " rad " << 180.0 / M_PI * phi_rad << " deg" << endl;
	phi_rad = M_PI / 2.0 - 0.5 * asin(tmp_arg);
	cout << "arcsin phi_1 " << phi_rad << " rad " << 180.0 / M_PI * phi_rad << " deg" << endl;
	phi_rad = atan((1.0 - sqrt(1.0 - tmp_arg*tmp_arg)) / tmp_arg);
	cout << "arctan phi_0 " << phi_rad << " rad " << 180.0 / M_PI * phi_rad << " deg" << endl;
	phi_rad = atan((1.0 + sqrt(1.0 - tmp_arg*tmp_arg)) / tmp_arg);
	cout << "arctan phi_1 " << phi_rad << " rad " << 180.0 / M_PI * phi_rad << " deg" << endl;
	#endif /* PROD */
	
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
	int const n_x = 1025;
	double const dx = (global_L - global_x_0) / (n_x - 1);
	double x;
	double Y[n_eq];
	Y[0] = global_y_0;  // Physical quantity y [m]
	Y[1] = global_v_0;  // Physical quantity v [m / s]
	Y[2] = phi;  // Physical quantity phi [rad]
	for (int i_x = 0 ; i_x < n_x; i_x++) {
		x = global_x_0 + i_x * dx;
		rungekutta4(x, dx, Y, rhs, n_eq);
	}
	return Y[0];  // Boundary value is y(L) = 0.
}
