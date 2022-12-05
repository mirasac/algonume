/*
MC things to do before PROD:
- Downgrade global variables which are not needed globally.
*/

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../../mclib/mclib.h"

static double const global_a = 0.25;
static double const global_k = 1.0;

double T_0(double x);

double f_L(double t);

double f_R(double t);

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	int const N_x = 256;  // Number of spatial points.
	int const N_t = 100;  // Number of temporal intervals.
	double const t_b = 0.0;
	double const t_e = 0.1;
	double const x_b = -1.0;
	double const x_e = 1.0;
	double const dt = (t_e - t_b) / N_t;
	double const dx = (x_e - x_b) / N_x;
	double const alpha = global_k * dt / (dx*dx);
	double const theta = 0.5;
	double x, t;
	double d_inf[N_x-2], d[N_x-2], d_sup[N_x-2];
	for (int i = 0; i <= N_x - 3; i++) {
		// MC finire con risoluzione tridiagonale.
	}
	// MC go ahead.
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
