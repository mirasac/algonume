#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

#define NUM_ITERATIONS 10
#define N_WIDTH 9

double function(double x) {
	return sin(x);
}

double function_1(double x) {
	return cos(x);
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	double x, h, epsilon_FD, epsilon_BD, epsilon_CD;
	ofstream plotfile;
	x = 1.0;
	h = 0.5;
	plotfile.open("plotfile.dat");
	plotfile.precision(N_PRECISION);
	//outfile.setf(ios::scientific);
	plotfile << setw(N_WIDTH) << "1/h " << "FD " << "BD " << "CD" << endl;
	for (int n = 0; n < NUM_ITERATIONS; n++) {
		epsilon_FD = fabs(forward_difference(function, h, x) - function_1(x));
		epsilon_BD = fabs(backward_difference(function, h, x) - function_1(x));
		epsilon_CD = fabs(central_difference(function, h, x) - function_1(x));
		plotfile << setw(N_WIDTH) << 1.0 / h << ' ' << epsilon_FD << ' ' << epsilon_BD << ' ' << epsilon_CD << endl;
		h /= 2;
	}
	plotfile.close();
	return 0;
}

