#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "../mclib.h"

#define NUM_INTERVALS 4

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
	// To avoid the singular point it is convenient sing the Taylor expansion up to an order of t where the machine epsilon is negligeble.
	if (t < 1e-3) {
		_return = 1.0 - t*t / (3.0 * 2) + t*t*t*t / 125.0;
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
	double m, t, x, x_0;
	double h[NUM_INTERVALS];
	x = 0.8;
	x_0 = 0.0;
	FILE * out = fopen("out.dat", "a");
	fputs("h, M(h), T(h)\n", out);
	for (int n = 0; n < NUM_INTERVALS; n++) {
		m = primitive(midpointquad, function, x, x_0, h[n]);
		t = primitive(trapezioidalquad, function, x, x_0, h[n]);
		fprintf(out, "%f, %f, %f\n", h[n], m, t);
	}
	return 0;
}

double primitive(double (*quad)(double (*f)(double x), double a, double b, int N), double (*f)(double x), double x, double x_0, double h) {
	int N;
	N = (x - x_0) / h; // We suppose that h is multiple of x - x_0.
	return quad(f, x_0, x, N);
}

