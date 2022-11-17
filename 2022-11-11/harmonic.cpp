#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

static double const global_omega = 2.0;
static double const global_mass = 1.0;

void a(double const X_0[], double R[]) {
	R[0] = - global_omega*global_omega * X_0[0];
}

void rhs(double const Y_0[], double R[]) {
	// MC continue.
}

double E(double x, double v) {
	// Mass is assumed 1.
	return 0.5 * global_mass * global_omega*global_omega * x*x + 0.5 * global_mass * v*v;
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	int const n_eq = 1;
	double const T = 2.0 * M_PI / global_omega;
	int N;
	double h, t, E_vv, E_rk2, t_min, t_max;
	double X_vv[n_eq], V_vv[n_eq], Y_rk2[n_eq];
	X_vv[0] = 1.0;
	V_vv[0] = 0.0;
	t_min = 0.0;
	t_max = 10 * T;
	E_vv = E(X_vv[0], V_vv[0]);
	Y_rk2[0] = X_vv[0];
	Y_rk2[1] = Y_vv[0];
	h = 0.02 * T;
	N = (t_max - t_min) / h;
	ofstream plot_file;
	plot_file.open("harmonic.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "t x_{vv}(t) v_{vv}(t) E_{vv} x_{rk2}(t) v_{rk2}(t) E_{rk2}" << endl;
	for (int i = 0; i <= N; i++) {
		t = i * h;
		plot_file << t;
		// Verlet results.
		plot_file << ' ' << X_vv[0] << ' ' << V_vv[0];
		plot_file << ' ' << fabs(E(X_vv[0], V_vv[0]) - E_vv) / E_vv;
		plot_file << ' ' << t / 8.0;  // MC cheat GNUplot: I need to scale in the right way the axis.
		verlet_velocity(h, X_vv, V_vv, a, n_eq);
		// Runge-Kutta results.
		plot_file << ' ' << Y_rk2[0] << ' ' << Y_rk2[1];
		plot_file << ' ' << fabs(E(Y_rk2[0], Y_rk2[1]) - E_rk2) / E_rk2 << endl;
		rungekutta2(t, Y_rk2, rhs, n_eq);
		plot_file << endl;
	}
	plot_file.close();
	return 0;
}

