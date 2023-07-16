/*
Pipeline to execute the program from terminal, CWD is root of algonume repo.
make MAIN=sistbiol/sistbiol
sistbiol/sistbiol.out
gnuplot
plot "sistbiol_solution.dat" using 2:3 index 0
*/

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

#define N_INIT 20
#define N_STEP 10

/* Global variables */
static double const global_alpha_1 = 1.0;
static double const global_alpha_2 = global_alpha_1;
static double const global_beta_1 = 200.0;
static double const global_beta_2 = 10.0;
static int const global_gamma_1 = 4;
static int const global_gamma_2 = global_gamma_1;
static double const global_K_1 = 30.0;
static double const global_K_2 = 1.0;
static double const global_v = 1.0;

void rhs(double const t, double const Y_0[], double R[]) {
	// Parameters.
	double v = Y_0[2];
	double omega_1 = pow(v * Y_0[1], global_gamma_1);
	double omega_2 = pow(Y_0[0], global_gamma_2);
	double omega_3 = pow(v * Y_0[4], global_gamma_1);
	double omega_4 = pow(Y_0[3], global_gamma_2);
	// Equations.
	R[0] = global_alpha_1 * (1.0 - Y_0[0]) - global_beta_1 * Y_0[0] * omega_1 / (global_K_1 + omega_1);
	R[1] = global_alpha_2 * (1.0 - Y_0[1]) - global_beta_2 * Y_0[1] * omega_2 / (global_K_2 + omega_2);
	R[2] = v;
	R[3] = global_alpha_1 * (1.0 - Y_0[3]) - global_beta_1 * omega_3 / (global_K_1 + omega_3);
	R[4] = global_alpha_2 * (1.0 - Y_0[4]) - global_beta_2 * omega_4 / (global_K_2 + omega_4);
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION) << scientific;
	int const n_eq = 5;
	double const x_min = 0.0;
	double const x_max = 1.0;
	double const y_min = 0.0;
	double const y_max = 1.0;
	double dx, x_0, dy, y_0;
	double t_min, t_max, dt, t;
	double Y[n_eq];
	ofstream plot_file;
	dx = (x_max - x_min) / N_INIT;
	dy = (y_max - y_min) / N_INIT;
	t_min = 0.0;
	t_max = 100.0;
	dt = (t_max - t_min) / N_STEP;
	plot_file.open("sistbiol_solution.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	Y[2] = global_v;  // Workaround to pass values to the RHS of the equation during the function call.
	for (int i_x = 0; i_x <= N_INIT; i_x++) {
		x_0 = i_x * dx;
		for (int i_y = 0; i_y <= N_INIT; i_y++) {
			y_0 = i_y * dy;
			Y[0] = x_0;
			Y[1] = y_0;
			Y[3] = x_0;
			Y[4] = y_0;
			plot_file << "t x(t) y(t) x_hill(t) y_hill(t)" << endl;
			for (int i_t = 0; i_t <= N_STEP; i_t++) {
				t = i_t * dt;
				plot_file << t << ' ' << Y[0] << ' ' << Y[1] << ' ' << Y[3] << ' ' << Y[4] << endl;
				eulerstep(t, dt, Y, rhs, n_eq);
			}
			plot_file << '\n' << endl;
		}
	}
	plot_file.close();
	return 0;
}
