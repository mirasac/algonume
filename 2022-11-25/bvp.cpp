#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

double rhs(double x) {
	return 1.0;
}

int main() {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	int const N = 30;
	double const y_lower = 0.0;  // Boundary condition y(0) = 0.
	double const y_upper = 0.0;  // Boundary condition y(1) = 0.
	double const x_min = 0.0;
	double const x_max = 1.0;
	double dx, x;
	double d_inf[N], d[N], d_sup[N], b[N], * y;
	// Assign values.
	dx = (x_max - x_min) / N;
	for (int i = 0; i < N; i++) {
		x = x_min + i * dx;
		d_inf[i] = 1.0;
		d[i] = -2.0;
		d_sup[i] = 1.0;
		b[i] = dx*dx * rhs(x);
	}
	b[0] -= y_lower;
	b[N-1] -= y_upper;
	// Open data file.
	ofstream plot_file;
	plot_file.open("bvp.dat");
	plot_file << "x y(x)" << endl;
	// Solve using tridiagonal matrix.
	y = tridiagonal_solver(d_inf, d, d_sup, b, N);
	plot_file << x << ' ' << y_lower;
	for (int i = 0; i < N; i++) {
		x = x_min + (i + 1) * dx;
		plot_file << x << ' ' << y[i] << endl;
	}
	plot_file << x << ' ' << y_upper;
	// Clean up.
	plot_file.close();
	delete[] y;
	return 0;
}

