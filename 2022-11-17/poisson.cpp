#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

void rhs(double r, double const Y_0[], double R[]) {
	R[0] = Y_0[1];
	R[1] = 0.5 * r * exp(-r);
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	ofstream plot_file;
	
	// Point 1.
	double const b = 20.0;
	double const s_step = 0.1;
	int const N = 1000;
	int const n_eq = 2;
	int const n_s = 0;  // MC just for test.
	double s, r_min, r_max, r, dr;
	double Y[n_eq];
	r_min = 0.0;
	r_max = b;
	dr = (r_max - r_min) / N;
	plot_file.open("poisson_1.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "r phi(r,s=0.0) phi(r,s=0.2) phi(r,s=0.4) phi(r,s=0.6) phi(r,s=0.8) phi(r,s=1.0)" << endl;
	// Integrate at different values of s.
	for (int i_s = 0; i_s <= n_s; i_s++) {
		s =  (2 * i_s) * s_step;
		Y[0] = 0.0;
		Y[1] = s;
		for (int i_r = 0; i_r <= N; i_r++) {
			r = r_min + i_r * dr;
			plot_file << r << ' ' << Y[0];
			rungekutta4(r, dr, Y, rhs, n_eq);
			plot_file << endl;
		}
	}
	plot_file.close();
	
	// Point 2.
	
	// Point 3.
	
	// Point 4.
	return 0;
}

