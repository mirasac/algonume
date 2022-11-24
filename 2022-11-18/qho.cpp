#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

static double const global_x_m = -1.0;
static double global_E = 0.5;

double solution(double x) {
	return exp(-0.5 * x*x);
}

double solution_1(double x) {
	return -x * solution(x);
}

double V(double x) {
	return 0.5 * x*x;
}

void rhs(double const x, double const Y_0[], double R[]) {
	R[0] = Y_0[1];
	R[1] = -2.0 * (global_E - V(x)) * Y_0[0];
}

double residual(double E) {
	global_E = E;
	int const N = 800;
	int const n_eq = 2;
	double x_a, x_b, x, dx, D, epsilon, Y_L, Y_L_1, Y_R, Y_R_1;
	double Y[n_eq];
	epsilon = 1e-10;
	x_a = -10.0;
	x_b = 10.0;
	D = 1.0;//sqrt(x_a*x_a + x_b*x_b) + epsilon;  // Commented his suggestion for the normalization.
	// Forward integration.
	dx = (global_x_m - x_a) / N;
	Y[0] = solution(x_a);
	Y[1] = solution_1(x_a);
	for (int n = 0; n <= N; n++) {
		x = x_a + n * dx;
		rungekutta4(x, dx, Y, rhs, n_eq);
	}
	Y_L = Y[0];
	Y_L_1 = Y[1];
	// Backward integration.
	dx = (global_x_m - x_b) / N;
	Y[0] = solution(x_b);
	Y[1] = solution_1(x_b);
	for (int n = 0; n <= N; n++) {
		x = x_b + n * dx;
		rungekutta4(x, dx, Y, rhs, n_eq);
	}
	Y_R = Y[0];
	Y_R_1 = Y[1];
	return (Y_L_1 * Y_R - Y_R_1 * Y_L) / D;
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	int const n_eq = 2;
	double x_a, x_b, x, dx;
	double Y[n_eq];
	ofstream plot_file;
	
	// Point 1.
	int const N = 800;
	plot_file.open("qho.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "x psi(x)" << endl;
	// Forward integration.
	x_a = -10.0;
	x_b = 10.0;
	dx = (x_b - x_a) / N;
	Y[0] = solution(x_a);
	Y[1] = solution_1(x_a);
	for (int n = 0; n <= N; n++) {
		x = x_a + n * dx;
		plot_file << x << ' ' << Y[0] << endl;
		rungekutta4(x, dx, Y, rhs, n_eq);
	}
	plot_file << '\n' << endl;
	// Backward integration.
	x_a = 10.0;
	x_b = -10.0;
	dx = (x_b - x_a) / N;
	Y[0] = solution(x_a);
	Y[1] = solution_1(x_a);
	for (int n = 0; n <= N; n++) {
		x = x_a + n * dx;
		plot_file << x << ' ' << Y[0] << endl;
		rungekutta4(x, dx, Y, rhs, n_eq);
	}
	plot_file.close();
	
	// Point 2.
	double E, E_min, E_max, dE, res;
	E_min = 0.0;
	E_max = 5.0;
	dE = (E_max - E_min) / N;
	plot_file.open("qho_residual.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "E residual(E)" << endl;
	for (int n = 0; n <= N; n++) {
		E = E_min + n * dE;
		res = residual(E);
		plot_file << E << ' ' << res << endl;
	}
	plot_file.close();
	
	// Point 3.
	cout << "\nPoint 3" << endl;
	int const n_intervals = 20;
	int n_bracketing;
	double const tollerance = 1e-6;
	double E_L[n_intervals], E_R[n_intervals];
	n_bracketing = bracketing(residual, E_min, E_max, n_intervals, E_L, E_R);
	for (int n = 0; n < n_bracketing; n++) {
		cout << bisection(residual, E_L[n], E_R[n], tollerance) << endl;
	}
	return 0;
}

