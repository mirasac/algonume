#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

/*
Formula of the bidimensional circle of radius r centered in (x_0, y_0): (x - x_0)^2 + (y - y_0)^2 = r^2
*/
double function(double x, double r, double x_0, double y_0) {
	return y_0 + sqrt(r*r - (x - x_0)*(x - x_0));
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	srand48(time(NULL));
	double tol, I, r, sigma, err, A_sq, x, y, f_x, f_y;
	int N, N_in;
	ofstream f_plot, f_err;
	r = 1.0;
	A_sq = 4.0 * r*r;
	tol = 1e-4;
	N = 4;
	f_plot.open("pi.dat");
	f_plot.precision(N_PRECISION);
	f_plot << "x " << "y" << endl;
	f_err.open("pi_err.dat");
	f_err.precision(N_PRECISION);
	f_err << "N " << "sigma" << endl;
	do {
		N_in = 0;
		sigma = 0.0;
		for (int i = 1; i <= N; ++i) {
			x = drand48() * 2.0 - 1.0;
			y = drand48() * 2.0 - 1.0;
			if (x*x + y*y <= r*r) {
				++N_in;
				f_plot << x << ' ' << y << endl;
				sigma += (M_PI - A_sq * N_in / N) * (M_PI - A_sq * N_in / N);
			}
		}
		// MC evaluation of sigma does not work.
		sigma /= (N - 1);
		f_err << N << ' ' << sigma << endl;
		I = A_sq * N_in / N;
		err = I / M_PI - 1.0;
		N *= 2;
	} while (fabs(err) > tol);
	f_plot.close();
	f_err.close();
	return 0;
}

