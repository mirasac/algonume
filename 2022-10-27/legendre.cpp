#include <iostream>
#include <iomanip>
#include "../mclib/mclib.h"

#define NUM_INTERVALS 100

double eval_legendre(double x) {
	int order;
	order = 13;
	return polynomial_legendre(order, x);
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	double a, b, tolerance;
	double x_L[NUM_INTERVALS], x_R[NUM_INTERVALS];
	int N;
	a = -1.0;
	b = 1.0;
	tolerance = 1e-15;
	N = bracketing(eval_legendre, a, b, NUM_INTERVALS, x_L, x_R);
	for (int n = 0; n < N; n++) {
		cout << "Root #" << n + 1 << ": " << bisection(eval_legendre, x_L[n], x_R[n], tolerance) << endl;
	}
	return 0;
}

