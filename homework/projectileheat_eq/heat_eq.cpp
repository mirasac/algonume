/*
MC things to do before PROD:
- Downgrade global variables which are not needed globally.
*/

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
#define DDD 0  // Plot in 3D, not part of the homework.

#define HOMEWORK_NAME "heat_eq"
#define N_PRECISION 8
#define TOLERANCE 1e-7

/* Global variables */
static double const global_a = 0.25;
static double const global_k = 1.0;

/* Function prototypes */
double T_0(double x);
double f_L(double t);
double f_R(double t);

int main() {
	// Set up.
	using namespace std;
	cout << setprecision(N_PRECISION) << scientific;
	int const N_x = 256;  // Number of spatial points.
	int const N_t = 100;  // Number of temporal intervals.
	double const t_e = 0.1;
	double const x_b = -1.0;
	double const x_e = 1.0;
	double const dt = t_e / N_t;
	double const dx = (x_e - x_b) / (N_x - 1);
	double const alpha = global_k * dt / (dx*dx);
	double const theta = 0.5;
	double x, t;
	double d_inf[N_x-2], d[N_x-2], d_sup[N_x-2], b[N_x-2], * T;
	ofstream plot_file;
	// Set tridiagonal matrix coefficients.
	for (int i_x = 0; i_x <= N_x-3; i_x++) {
		d_inf[i_x] = - alpha * theta;
		d[i_x] = 2.0 * alpha * theta + 1.0;
		d_sup[i_x] = - alpha * theta;
	}
	// Set initial values.
	T = new double[N_x-2];
	for (int i_x = 0; i_x <= N_x - 3; i_x++) {
		x = x_b + (i_x + 1) * dx;
		T[i_x] = T_0(x);
	}
	plot_file.open(HOMEWORK_NAME ".dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "x T(x)" << endl;
	// Loop over time.
	for (int i_t = 1; i_t <= N_t; i_t++) {
		t = i_t * dt;
		// Print solutions at requested time steps.
		if (i_t == N_t / 2 || i_t == N_t || DDD) {
			plot_file << x_b;
			#if DDD == 1
			plot_file << ' ' << t;
			#endif /* DDD == 1 */
			plot_file << ' ' << f_L(0.0) << endl;
			for (int i = 0; i <= N_x - 3; i++) {
				x = x_b + (i + 1) * dx;
				plot_file << x;
				#if DDD == 1
				plot_file << ' ' << t;
				#endif /* DDD == 1 */
				plot_file << ' ' << T[i] << endl;
			}
			plot_file << x_e;
			#if DDD == 1
			plot_file << ' ' << t;
			#endif /* DDD == 1 */
			plot_file << ' ' << f_R(0.0) << endl;
			#if DDD == 1
			plot_file << endl;
			#else
			plot_file << '\n' << endl;
			// MC fare meglio il grafico, magari separarlo su stesso file ma due colonne dati diverse.
			#endif /* DDD == 1 */
			cout << "Printed solutions at t = " << t << endl;
		}
		// Loop over space at fixed time step.
		b[0] = alpha * (1.0 - theta) * f_L(0.0) + (1.0 - 2.0 * alpha * (1.0 - theta)) * T[0] + alpha * (1.0 - theta) * T[1];
		for (int i_x = 1; i_x <= N_x - 4; i_x++) {
			b[i_x] = alpha * (1.0 - theta) * T[i_x-1] + (1.0 - 2.0 * alpha * (1.0 - theta)) * T[i_x] + alpha * (1.0 - theta) * T[i_x+1];
		}
		b[N_x-3] = alpha * (1.0 - theta) * T[N_x-4] + (1.0 - 2.0 * alpha * (1.0 - theta)) * T[N_x-3] + alpha * (1.0 - theta) * f_R(0.0);
		delete[] T;
		T = tridiagonal_solver(d_inf, d, d_sup, b, N_x-2);
	}
	// Plot values at desired times.
	// Tear down.
	delete[] T;
	plot_file.close();
	return 0;
}

double T_0(double x) {
	return 1.0 + exp(-x*x / global_a*global_a);
}

double f_L(double t) {
	return 1.0;
}

double f_R(double t) {
	return 1.0;
}
