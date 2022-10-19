#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <iomanip>
#include "../mclib/mclib.h"

#define NUM_INTERVALS 4
#define PRECISION 8
#define COL_WIDTH 10

/*
Mathematical function I can change without editing the code for the quadrature.
*/
double function(double x) {
	double _return;
	/*
	if (x != 0.0) {
		_return = sin(x) / x;
	} else {
		_return = 1.0;
	}
	*/
	// To avoid the singular point it is convenient sing the Taylor expansion up to an order of t where the machine epsilon is negligible.
	if (x < 1e-3) {
		_return = 1.0 - x*x / (3.0 * 2) + x*x*x*x / 125.0;
	} else {
		_return = sin(x) / x;
	}
	return _return;
}

double primitive(
	double (*quad)(
		double (*f)(double x),
		double a,
		double b,
		int N
	),
	double (*f)(double x),
	double x,
	double x_0,
	double h
);

int main() {
	using namespace std;
	double m, t, x, x_0, h, a, b, delta_x, delta_h;
	double h_interval[NUM_INTERVALS], sum_prev[NUM_INTERVALS];
	int num_x;
	
	// Datafile with integrals at various interval divisions.
	x = 0.8;
	x_0 = 0.0;
	ofstream outfile;
	outfile.open("outfile.dat");
	outfile.precision(PRECISION);
	//outfile.setf(ios::scientific);
	outfile << setw(COL_WIDTH) << "h " << "M(h) " << "T(h)" << endl;
	for (int n = 0; n < NUM_INTERVALS; n++) {
		h = (x - x_0) / (1 << n);
		m = primitive(midpointquad, function, x, x_0, h);
		t = primitive(trapezioidalquad, function, x, x_0, h);
		outfile << setw(COL_WIDTH) << h << " " << m << " " << t << endl;
	}
	outfile.close();
	
	// Datafile for plot of sin(x) / x with various interval divisions.
	// Input data.
	a = 0.0;
	b = 25.0;
	num_x = 1000;
	// Output file preparation.
	ofstream plotfile;
	plotfile.open("plotfile.dat");
	plotfile.precision(PRECISION);
	plotfile << setw(COL_WIDTH) << "#x";
	// Print number of intervals in the header.
	for (int n = 0; n < NUM_INTERVALS; n++) {
		plotfile << setw(COL_WIDTH) << " " << n+1;
		sum_prev[n] = 0.0;
	}
	plotfile << endl;
	// Integral evaluation.
	orderinterval(&a, &b);
	delta_x = (b - a) / num_x;
	x_0 = a;
	for (int i = 1; i <= num_x; i++) {
		x = a + delta_x * i;
		delta_h = x - x_0;
		// MC mettere setw(COL_WIDTH) dove serve.
		plotfile << x;
		for (int n = 0; n < NUM_INTERVALS; n++) {
			h = delta_h / (1 << n);
			sum_prev[n] += primitive(trapezioidalquad, function, x, x_0, h);
			plotfile << " " << sum_prev[n];
		}
		plotfile << endl;
		x_0 = x;
	}
	outfile.close();
	return 0;
}

double primitive(double (*quad)(double (*f)(double x), double a, double b, int N), double (*f)(double x), double x, double x_0, double h) {
	int N;
	N = (x - x_0) / h; // We suppose that h is multiple of x - x_0.
	return quad(f, x_0, x, N);
}

