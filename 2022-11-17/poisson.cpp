#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

void rhs(double r, double const Y_0[], double R[]) {
	R[0] = Y_0[1];
	R[1] = -0.5 * r * exp(-r);
}

double solution(double r) {
	return 1.0 - 0.5 * (r + 2) * exp(-r);
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
	double r_min, r_max, r, dr;
	double Y_1[n_eq], Y_2[n_eq], Y_3[n_eq], Y_4[n_eq], Y_5[n_eq], Y_6[n_eq];
	r_min = 0.0;
	r_max = b;
	dr = (r_max - r_min) / N;
	Y_1[0] = 0.0;
	Y_1[1] = 0.0;
	Y_2[0] = 0.0;
	Y_2[1] = 0.2;
	Y_3[0] = 0.0;
	Y_3[1] = 0.4;
	Y_4[0] = 0.0;
	Y_4[1] = 0.6;
	Y_5[0] = 0.0;
	Y_5[1] = 0.8;
	Y_6[0] = 0.0;
	Y_6[1] = 1.0;
	plot_file.open("poisson_1.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "r phi(r,s=0.0) phi(r,s=0.2) phi(r,s=0.4) phi(r,s=0.6) phi(r,s=0.8) phi(r,s=1.0)" << endl;
	for (int i_r = 0; i_r <= N; i_r++) {
		r = r_min + i_r * dr;
		plot_file << r;
		plot_file << ' ' << Y_1[0];
		rungekutta4(r, dr, Y_1, rhs, n_eq);
		plot_file << ' ' << Y_2[0];
		rungekutta4(r, dr, Y_2, rhs, n_eq);
		plot_file << ' ' << Y_3[0];
		rungekutta4(r, dr, Y_3, rhs, n_eq);
		plot_file << ' ' << Y_4[0];
		rungekutta4(r, dr, Y_4, rhs, n_eq);
		plot_file << ' ' << Y_5[0];
		rungekutta4(r, dr, Y_5, rhs, n_eq);
		plot_file << ' ' << Y_6[0];
		rungekutta4(r, dr, Y_6, rhs, n_eq);
		plot_file << endl;
	}
	plot_file.close();
	
	// Point 2.
	
	// Point 3.
	
	// Point 4.
	plot_file.open("poisson_4.dat");
	plot_file << setprecision(N_PRECISION) << scientific;
	plot_file << "r phi(r)" << endl;
	for (int i_r = 0; i_r <= N; i_r++) {
		r = r_min + i_r * dr;
		plot_file << r << ' ' << solution(r) << endl;
	}
	plot_file.close();
	return 0;
}

