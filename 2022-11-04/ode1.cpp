// set key autotitle columnhead

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

#define N_STEP_SIZE 9

double solution(double x) {
	return exp(-x*x / 2.0);
}

void rhs(double const t, double const Y_0[], double R[]) {
	R[0] = -t * Y_0[0];
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION) << scientific;
	int N, h_int;
	int const n_eq = 1;
	double t_min, t_max, h, t, y_ref, abs_err, rel_err, h_order;
	double Y[n_eq];
	char * filename = new char[9 + N_STEP_SIZE / 10 + 1];
	ofstream plot_file;
	t_min = 0.0;
	t_max = 3.0;
	h_int = 5;
	h_order = 0.1;
	for (int n = 1; n <= N_STEP_SIZE; ++n) {
		h = h_int * h_order;
		sprintf(filename, "ode1_%d.dat", n);
		plot_file.open(filename);
		plot_file << setprecision(N_PRECISION) << scientific;
		plot_file << "t y(t) abs_{err}(t) rel_{err}(t)" << endl;
		Y[0] = 1.0;  // Initial value.
		N = (t_max - t_min) / h;  // Suppose h multiple of t_max - t_min.
		for (int i = 0; i <= N; ++i) {
			t = i * h;
			y_ref = solution(t);
			abs_err = fabs(Y[0] - y_ref);
			rel_err = abs_err / fabs(y_ref);
			plot_file << t << ' ' << Y[0] << ' ' << abs_err << ' ' << rel_err << endl;
			eulerstep(t, h, Y, rhs, n_eq);
		}
		plot_file.close();
		h_int /= 2;
		if (n % 3 == 0) {
			h_order /= 10;
			h_int = 5;
		}
	}
	delete[] filename;
	return 0;
}
