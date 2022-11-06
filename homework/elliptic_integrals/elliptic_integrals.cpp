#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#ifndef PROD
#include <time.h>
#include "../../mclib/mclib.h"
#endif /* PROD */

#ifndef N_PRECISION
#define N_PRECISION 12
#endif /* N_PRECISION */

#define TOLLERANCE 1e-7
#define HOMEWORK_NAME "elliptic_integrals"

/* Global variables */
static double const global_k = 0.8;
static double const global_phi_c = M_PI / 4.0;
static double global_t = 0.0;

/* Functions */
double g(double u) {
	return 1.0 / sqrt(1.0 - global_k*global_k * sin(u)*sin(u));
}

#ifndef PROD
double g_approx(double u) {
	return 1.0 / sqrt(1.0 - global_k*global_k * u*u);
}
#endif /* PROD */

// This function is valid if -1/k <= phi <= 1/k.
double F_small_angles(double phi) {
	return asin(global_k * phi) / global_k;
}

double F(double phi) {
	double _return, phi_0;
	phi_0 = 0.0;
	if (fabs(phi - phi_0) <= TOLLERANCE) {
		_return = 0.0;
	} else {
		_return = gaussquad(g, phi_0, phi, 100, 4);
	}
	return _return;
}

double g1(double u) {
	double g_temp = g(u);
	return global_k*global_k / 2.0 * sin(2.0 * u) * g_temp*g_temp*g_temp;
}

double g2(double u) {
	double g_temp = g(u);
	double cos_temp = cos(2.0 * u);
	return global_k*global_k * g_temp*g_temp*g_temp * ( cos_temp + 3.0 / 4.0 * global_k*global_k * (1.0 - cos_temp*cos_temp) * g_temp*g_temp );
}

double h_0(double phi) {
	return F(global_phi_c) - global_t;
}

double h_1(double phi) {
	return h_0(phi) + (phi - global_phi_c) * g(global_phi_c);
}

double h_2(double phi) {
	double delta_phi = phi - global_phi_c;
	return h_1(phi) + delta_phi * g1(global_phi_c) / 2.0 * delta_phi;
}

double h_3(double phi) {
	double delta_phi = phi - global_phi_c;
	return h_2(phi) + delta_phi * g2(global_phi_c) / 6.0 * delta_phi*delta_phi;
}

/* Entry point */
int main() {
	// Setup.
	using namespace std;
	cout << setprecision(N_PRECISION);
	
	// Point 1.
	cout << "Point 1" << endl;
	double T, a, b;
	int N, Ng;
	N = 100;
	Ng = 4;
	a = 0.0;
	b = M_PI / 2.0;
	// Assume L / g = 1.
	T = 4.0 * rectangularquad(g, a, b, N);
	cout << "Rectangular: T = " << T << endl;
	T = 4.0 * midpointquad(g, a, b, N);
	cout << "Midpoint: T = " << T << endl;
	T = 4.0 * trapezioidalquad(g, a, b, N);
	cout << "Trapezioidal: T = " << T << endl;
	T = 4.0 * simpsonquad(g, a, b, N);
	cout << "Simpson: T = " << T << endl;
	T = 4.0 * gaussquad(g, a, b, N, Ng);
	cout << "Gauss: T = " << T << endl;

	// Point 2.
	cout << "\nPoint 2" << endl;
	double T_small_angles;
	// Assume L / g = 1.
	T_small_angles = 2.0 * M_PI;
	cout << "T_small_angles = " << T_small_angles << endl;
	cout << "Relative error of small angle approximation with respect to real period (evaluated with Gauss method): " << fabs(T - T_small_angles) / T << endl;

	// Point 3.
	cout << "\nPoint 3" << endl;
	#ifndef PROD
	clock_t time_start, time_end;
	time_start = clock();
	#endif /* PROD */
	double phi_t;
	ofstream file_plot;
	file_plot.open(HOMEWORK_NAME ".dat");
	file_plot << setprecision(N_PRECISION);
	file_plot << "t" << " phi(t)";
	#ifndef PROD
	file_plot << " phi_1(t)" << " phi_2(t)" << " phi_3(t)" << endl;
	#endif /* PROD */
	file_plot << 0.0 << ' ' << 0.0; // Results evaluated analytically.
	#ifndef PROD
	file_plot << ' ' << 0.0 << ' ' << 0.0 << ' ' << 0.0; // Results evaluated analytically.
	#endif /* PROD */
	file_plot << endl;
	b = 100.0; // Conservative value to be able to invert the function to the maximum t requested.
	// Loop on t;
	for (int i = 1; i <= 300; ++i) {
		global_t = i * 0.1;
		file_plot << global_t;
		phi_t = bisection(F, a, b, TOLLERANCE);
		file_plot << ' ' << phi_t;
		#ifndef PROD
		phi_t = bisection(h_1, a, b, TOLLERANCE);
		file_plot << ' '  << phi_t;
		phi_t = bisection(h_2, a, b, TOLLERANCE);
		file_plot << ' '  << phi_t;
		phi_t = bisection(h_3, a, b, TOLLERANCE);
		file_plot << ' '  << phi_t;
		#endif /* PROD */
		file_plot << endl;
		// MC fare calcoli anche con Newton.
	}
	#ifndef PROD
	time_end = clock();
	cout << "CPU time employed for the operations: " << ((double) (time_end - time_start)) / CLOCKS_PER_SEC << " s" << endl;
	#endif /* PROD */

	// Exit.
	file_plot.close();
	return 0;
}
