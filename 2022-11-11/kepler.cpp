#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

void rhs(double t, double Y_0[], double R[]) {
	r = sqrt(Y_0[1]*Y_0[1] + Y_0[3]*Y_0[3]);
	// Assume GM = 1.
	R[0] = Y_0[1];
	R[1] = - Y_0[1] / (r*r*r);
	R[2] = Y_0[3];
	R[3] = - Y_0[3] / (r*r*r);
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	// MC go ahead.
	return 0;
}

