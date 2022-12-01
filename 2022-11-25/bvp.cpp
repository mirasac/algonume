#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

int main() {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	int const N = 32;  // Points, two of which are boundaries.
	double const x_min = 0.0;
	double const x_max = 1.0;
	double const y_lower = 0.0;  // Boundary condition y(0) = 0.
	double const y_upper = 0.0;  // Boundary condition y(1) = 0.
	double dx, x;
	double d_inf[N], d[N], d_sup[N], b[N], * y;
	dx = (x_max - x_min) / (N - 1);
	// Initialize arrays.
	for (int i = 0; i < N; i++) {
		d_inf[i] = 1.0;
		d[i] = -2.0;
		d_sup[i] = 1.0;
		b[i] = dx*dx;  // RHS is set to 1.
	}
	b[1] -= y_lower;
	b[N-2] -= y_upper;
	// Open data file.
	ofstream plot_file;
	plot_file.open("bvp.dat");
	plot_file << "x y(x)" << endl;
	// Solve using tridiagonal matrix.
	y = tridiagonal_solver(d_inf+1, d+1, d_sup+1, b+1, N-2);
	cout << y[0] << ' ' << y[N-1] << endl;
	plot_file << x_min << ' ' << y_lower << endl;
	for (int i = 0; i <= N-3; i++) {
		x = x_min + (i + 1) * dx;
		plot_file << x << ' ' << y[i] << endl;
	}
	plot_file << x_max << ' ' << y_upper << endl;
	// Clean up.
	plot_file.close();
	delete[] y;
	return 0;
}
