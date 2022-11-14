#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

static double const global_omega = 6.0;

void a(double const X_0[], double R[]) {
	R[0] = - global_omega*global_omega * X_0[0];
}

double E(double x, double v) {
	// Mass is assumed 1.
	return 0.5 * global_omega*global_omega * x*x + 0.5 * v*v;
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	int const n_eq = 1;
	double const T = 2.0 * M_PI / global_omega;
	int N;
	double h, t, E_vv, E_rk2, t_min, t_max;
	double X_vv[n_eq], V_vv[n_eq], X_rk2[n_eq], V_rk2[n_eq];
	X_vv[0] = 1.0;
	V_vv[0] = 0.0;
	t_min = 0.0;
	t_max = 10 * T;
	E_vv = E(X_vv[0], V_vv[0]);
	//Y_rk2[0] = 1.0;
	//V_rk2[0] = 0.0;
	h = 0.02 * T;
	N = (t_max - t_min) / h;
	ofstream plot_file;
	plot_file.open("harmonic.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	for (int i = 0; i <= N; i++) {
		t = i * h;
		plot_file << t << ' ' << fabs(E(X_vv[0], V_vv[0]) - E_vv) / E_vv << endl;
		verlet_velocity(h, X_vv, V_vv, a, n_eq);
	}
	plot_file.close();
	return 0;
}

