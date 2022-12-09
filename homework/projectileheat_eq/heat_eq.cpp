#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

// Do not copy into report.
#include "../../mclib/mclib.h"
#ifdef N_PRECISION
#undef N_PRECISION
#endif /* N_PRECISION */

#define HOMEWORK_NAME "heat_eq"
#define N_PRECISION 8

/* Function prototypes */
double T_0(double x);
double f_L(double t);
double f_R(double t);

int main() {
	// Set up.
	using namespace std;
	cout << setprecision(N_PRECISION) << scientific;
	// Parameters.
	int const N_x = 256;
	int const N_t = 100;
	double const t_e = 0.1;
	double const x_b = -1.0;
	double const x_e = 1.0;
	double const k = 1.0;
	double const theta = 0.5;
	double const dt = t_e / N_t;
	double const dx = (x_e - x_b) / (N_x - 1);
	double const alpha = k * dt / (dx*dx);
	double const c_1 = - alpha * theta;
	double const c_2 = 1.0 - 2.0 * c_1;
	// Variables.
	double t;
	double x[N_x], d_inf[N_x], d[N_x], d_sup[N_x], b[N_x], T[N_x];
	ofstream plot_file;

	// Set spatial grid and tridiagonal matrix coefficients.
	for (int i = 0; i < N_x; i++) {
		x[i] = x_b + i * dx;
		d_inf[i] = c_1;
		d[i] = c_2;
		d_sup[i] = c_1;
	}

	// Assign initial and boundary values.
	T[0] = f_L(0.0);
	for (int i = 1; i <= N_x - 2; i++) {
		T[i] = T_0(x[i]);
	}
	T[N_x-1] = f_R(0.0);

	// Loop over time to evolve temperature distribution.
	plot_file.open(HOMEWORK_NAME ".dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	for (int n = 0; n <= N_t; n++) {
		t = n * dt;
		// Save solutions to file.
		if (n == N_t / 2 || n == N_t) {
			for (int i = 0; i < N_x; i++) {
				plot_file << x[i] << ' ' << T[i] << endl;
			}
			plot_file << '\n' << endl;
			cout << "Temperature distribution at t = " << t << " saved to file" << endl;
		}
		// Assign values to constant terms of the linear system.
		for (int i = 1; i <= N_x - 2; i++) {
			b[i] = (alpha + c_1) * T[i-1] + (c_2 - 2.0 * alpha) * T[i] + (alpha + c_1) * T[i+1];
		}
		b[1] -= d_inf[0] * f_L(t + dt);
		b[N_x-2] -= d_sup[N_x-1] * f_R(t + dt);
		// Evaluate T(x) at next time step.
		tridiagonal_solver_2(d_inf+1, d+1, d_sup+1, b+1, T+1, N_x-2);
	}

	// Teardown.
	plot_file.close();
	return 0;
}

double T_0(double x) {
	return 1.0 + exp(-x*x * 16.0);
}

double f_L(double t) {
	return 1.0;
}

double f_R(double t) {
	return 1.0;
}
