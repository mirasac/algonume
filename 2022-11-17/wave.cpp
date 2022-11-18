#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

static double global_k = 1.0;

void rhs(double const t, double const Y_0[], double R[]) {
	R[0] = Y_0[1];
	R[1] = - global_k*global_k * Y_0[0];
}

double residual(double k) {
	int const N = 1000;
	int const n_eq = 2;
	double x, dx, x_min, x_max;
	double Y[2];
	x_min = 0.0;
	x_max = 1.0;
	dx = (x_max - x_min) / N;
	Y[0] = 0.0;
	Y[1] = 1.0;  // This is s.
	global_k = k;  // MC probabilmente non funziona, passarlo come terzo valore a Y.
	for (int n = 0; n <= N; n++) {
		x = x_min + n * dx;
		rungekutta4(x, dx, Y, rhs, n_eq);
	}
	return Y[0];  // Final boundary condition is phi(1) = 0.
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	int const N = 100;
	int const n_eq = 2;
	double x_min, x_max, x, dx, s;
	double Y[n_eq];
	ofstream plot_file;
	x_min = 0.0;
	x_max = 1.0;
	dx = (x_max - x_min) / N;
	s = 1.0;
	
	// Point 1.
	global_k = 1.0;
	Y[0] = 0.0;
	Y[1] = s;
	plot_file.open("wave_1.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "x phi(x)" << endl;
	for (int n = 0; n <= N; n++) {
		x = x_min + n * dx;
		plot_file << x << ' ' << Y[0] << endl;
		rungekutta4(x, dx, Y, rhs, n_eq);
	}
	plot_file.close();
	
	// Point 2.
	plot_file.open("wave_2.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "x phi(x)" << endl;
	for (int i = 1; i <= 5; i++) {
		global_k = i;
		Y[0] = 0.0;
		Y[1] = s;
		for (int n = 0; n <= N; n++) {
			x = x_min + n * dx;
			plot_file << x << ' ' << Y[0] << endl;
			rungekutta4(x, dx, Y, rhs, n_eq);
		}
		plot_file << '\n' << endl;
	}
	plot_file.close();
	
	// Point 3.
	cout << "\nPoint 3" << endl;
	double tollerance;
	tollerance = 1e-6;
	x_min = 3.0;
	x_max = 4.0;
	cout << bisection(residual, x_min, x_max, tollerance) << endl;
	
	// Point 4.
	cout << "\nPoint 4" << endl;
	x_min = 0.0;
	x_max = 20.0;
	int const n_brack_max = 8;
	int n_brack;
	double x_L[n_brack], x_R[n_brack];
	n_brack = bracketing(residual, x_min, x_max, n_brack_max, x_L, x_R);
	for (int n = 0; n < n_brack; n++) {
		cout << bisection(residual, x_L[n], x_R[n], tollerance) << endl;
	}
	return 0;
}

