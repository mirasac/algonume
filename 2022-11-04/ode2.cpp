#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

#define N_STEP 200

void solution(double x, double Y[]) {
	// MC evaluate analytic solutions.
	Y[0] = cos(x);
	Y[1] = sin(x);
}

void rhs(double t, double Y_0[], double R[]) {
	R[0] = Y_0[1];
	R[1] = -Y_0[0];
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION) << scientific;
	
	// First part.
	int const n_eq = 2;
	double t_min, t_max, h, t, y_ref;
	double Y_ref[n_eq], Y_em[n_eq], Y_rk2[n_eq],  Y_rk4[n_eq];
	//char * filename = new char[9 + N_STEP_SIZE / 10 + 1];
	ofstream plot_file;
	t_min = 0.0;
	t_max = 2.0 * M_PI;
	Y_em[0] = 1.0;
	Y_em[1] = 0.0;
	Y_rk2[0] = 1.0;
	Y_rk2[1] = 0.0;
	Y_rk4[0] = 1.0;
	Y_rk4[1] = 0.0;
	h = (t_max - t_min) / N_STEP;
	plot_file.open("ode2.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "t x(t) y(t) x_em(t) y_em(t) x_rk2(t) y_rk2(t) x_rk4(t) y_rk4(t)" << endl;
	for (int i = 0; i <= N_STEP; i++) {
		//sprintf(filename, "ode1_%d.dat", n);
		t = i * h;
		eulerstep(t, h, Y_em, rhs, n_eq);
		rungekutta2(t, h, Y_rk2, rhs, n_eq);
		rungekutta4(t, h, Y_rk4, rhs, n_eq);
		solution(t, Y_ref);
		plot_file << t << ' ' << Y_ref[0] << ' ' << Y_ref[1] << ' ' << Y_em[0] << ' ' << Y_em[1] << ' ' << Y_rk2[0] << ' ' << Y_rk2[1] << ' ' << Y_rk4[0] << ' ' << Y_rk4[1] << endl;
	}
	plot_file.close();
	
	// Second part.
	double err;
	for (int N = 4; N <= 2048; N *= 2) {
		err = fabs();  // MC finire.
	} 
	//delete[] filename;
	return 0;
}
