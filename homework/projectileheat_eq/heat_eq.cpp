/*
MC things to do before PROD:
- Downgrade global variables which are not needed globally.
- use as indexes n for time and i for space.
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
#define TOLERANCE 1e-7
#define N_PRECISION 8

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
	double d_inf[N_x], d[N_x], d_sup[N_x], b[N_x], T[N_x], * T_tmp;
	ofstream plot_file;
	// Set tridiagonal matrix coefficients.
	for (int i_x = 0; i_x < N_x; i_x++) {
		d_inf[i_x] = - alpha * theta;
		d[i_x] = 2.0 * alpha * theta + 1.0;
		d_sup[i_x] = - alpha * theta;
	}
	// Assign initial and boundary values.
	T_tmp = new double[N_x-2];  // Needed as support array because of the implementation of tridiagonal_solver.
	T[0] = f_L(0.0);
	for (int i_x = 1; i_x <= N_x - 2; i_x++) {
		x = x_b + i_x * dx;
		T[i_x] = T_0(x);
		T_tmp[i_x-1] = T[i_x];
	}
	T[N_x-1] = f_R(0.0);
	plot_file.open(HOMEWORK_NAME ".dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "x T(x)" << endl;
	// Loop over time.
	for (int i_t = 0; i_t <= N_t; i_t++) {
		t = i_t * dt;
		// Print solutions at requested time steps.
		if (i_t == N_t / 2 || i_t == N_t || DDD || i_t == 0 || 1) {
			for (int i = 0; i < N_x; i++) {
				x = x_b + i * dx;
				plot_file << x;
				#if DDD == 1
				plot_file << ' ' << t;
				#endif /* DDD == 1 */
				plot_file << ' ' << T[i] << endl;
			}
			#if DDD == 1
			plot_file << endl;
			#else
			plot_file << '\n' << endl;
			// MC fare meglio il grafico, magari separarlo su stesso file ma due colonne dati diverse.
			#endif /* DDD == 1 */
			cout << "Printed solutions at t = " << t << endl;
		}
		// Assign values to constant terms of the linear system.
		//b[0] = alpha * (1.0 - theta) * f_L(0.0) + (1.0 - 2.0 * alpha * (1.0 - theta)) * T[0] + alpha * (1.0 - theta) * T[1];  // MC old version.
		for (int i_x = 1; i_x <= N_x - 2; i_x++) {
			T[i_x] = T_tmp[i_x-1];
			b[i_x] = alpha * (1.0 - theta) * T[i_x-1] + (1.0 - 2.0 * alpha * (1.0 - theta)) * T[i_x] + alpha * (1.0 - theta) * T[i_x+1];
		}
		//b[N_x-1] = alpha * (1.0 - theta) * T[N_x-2] + (1.0 - 2.0 * alpha * (1.0 - theta)) * T[N_x-1] + alpha * (1.0 - theta) * f_R(0.0);  // MC old version.
		b[1] -= d_inf[0] * f_L(t + dt);
		b[N_x-2] -= d_sup[N_x-1] * f_R(t + dt);
		delete[] T_tmp;
		// Evaluate T(x) at next time step.
		T_tmp = tridiagonal_solver(d_inf+1, d+1, d_sup+1, b+1, N_x-2);
	}
	// Tear down.
	delete[] T_tmp;
	plot_file.close();
	return 0;
}

double T_0(double x) {
	return 1.0 + exp(-x*x / (global_a*global_a));
}

double f_L(double t) {
	return 1.0;
}

double f_R(double t) {
	return 1.0;
}
