#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

double function(double x, double mean, double sigma) {
	return 1.0 / (sigma * sqrt(2.0 * M_PI)) * exp(-0.5 * ((x - mean) / sigma)*((x - mean) / sigma));
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	srand48(time(NULL));
	double mean, sigma, x, y, f, C;
	int N;
	C = 1.0; // In a general case should be a function of the point at which f is evaluated.
	N = 100000;
	mean = 0.0;
	sigma = 0.5;
	ofstream f_plot;
	f_plot.open("gaussian.dat");
	f_plot.precision(N_PRECISION);
	f_plot << "x " << "y" << endl;
	for (int k = 0; k < N; k++) {
		x = drand48() * 10.0 - 5.0;
		y = drand48() * C;
		f = function(x, mean, sigma);
		if (y <= f) {
			f_plot << x << ' ' << y << endl;
		}
	}
	f_plot.close();
	return 0;
}

