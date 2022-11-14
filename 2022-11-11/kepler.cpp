#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

static double const global_alpha = 0.3;  // Bounded orbits if < 2.

void rhs(double const t, double const Y_0[], double R[]) {
	double const r = sqrt(Y_0[0]*Y_0[0] + Y_0[2]*Y_0[2]);
	// Assume GM = 1.
	R[0] = Y_0[1];
	R[1] = - Y_0[0] / (r*r*r);
	R[2] = Y_0[3];
	R[3] = - Y_0[2] / (r*r*r);
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	int const n_eq = 4;
	int const N = 200;
	double t_min, t_max, t, dt;
	double Y[n_eq];
	ofstream plot_file;
	t_min = 0.0;
	t_max = 40.0;
	t = t_min;
	dt = (t_max - t_min) / N;
	Y[0] = 4.0;
	Y[1] = 0.0; // Initial velocity in y direction.
	Y[2] = 0.0;
	Y[3] = sqrt(global_alpha / sqrt(Y[0]*Y[0] + Y[2]*Y[2]));
	plot_file.open("kepler.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "t x(t) y(t)" << endl;
	for (int i = 1; i <= N; i++) {
		rungekutta4(t, dt, Y, rhs, n_eq);
		t = i * dt;
		//cout << t << ' ' << Y[0] << ' ' << Y[2] << endl;
		plot_file << t << ' ' << Y[0] << ' ' << Y[2] << endl;
	}
	plot_file.close();
	return 0;
}

