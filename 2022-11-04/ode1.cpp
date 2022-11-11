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

void rhs(double t, double Y_0[], double R[]) {
	R[0] = -t * Y_0[0];
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION) << scientific;
	int N;
	int const n_eq = 1;
	double t_min, t_max, h, t, y_ref, abs_err, rel_err;
	double Y[n_eq];
	char * filename = new char[9 + N_STEP_SIZE / 10 + 1];
	ofstream plot_file;
	t_min = 0.0;
	t_max = 3.0;
	h = 0.5;
	for (int n = 1; n <= N_STEP_SIZE; ++n) {
		sprintf(filename, "ode1_%d.dat", n);
		plot_file.open(filename);
		plot_file << setprecision(N_PRECISION) << scientific;
		plot_file << "t y(t) abs_err rel_err" << endl;
		Y[0] = 1.0;  // Initial value.
		N = (t_max - t_min) / h;  // Suppose h multiple of t_max - t_min.
		for (int i = 0; i <= N; ++i) {
			t = i * h;
			eulerstep(t, h, Y, rhs, n_eq);
			y_ref = solution(t);
			abs_err = fabs(Y[0] - y_ref);
			rel_err = abs_err / fabs(y_ref);
			plot_file << t << ' ' << Y[0] << ' ' << abs_err << ' ' << rel_err << endl;
		}
		plot_file.close();
		h /= 2;  // MC finire perché non step corretti.
		cout << h << endl; // MC debug.
	}
	delete[] filename;
	return 0;
}
