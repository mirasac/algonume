#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

double solution(double x) {
	return exp(-x*x / 2.0);
}

double rhs(double t, double Y_0[], double R[]) {
	R[0] = -t * Y_0[0];
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	double y_0, y_min, y_max;
	ofstream plot_file;
	y_0 = 1.0;
	y_min = 0.0;
	y_max = 3.0;
	plot_file.open("ode1.dat");
	plot_file << setprecision(N_PRECISION);
	// MC go ahead, page 9 chapter 7.
	plot_file.close();
	return 0;
}
