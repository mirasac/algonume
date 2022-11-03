#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

#define NUM_ITERATIONS 1000
#define MAX_ITERATIONS 1000000

int main() {
	// Setup.
	using namespace std;
	srand48(time(NULL));
	cout << setprecision(N_PRECISION);
	// Variables declaration.
	ofstream f_distribution, f_moments;
	double k_1, k_2, x_i;
	int N;
	// Loop and file write for distribution plot.
	f_distribution.open("distribution.dat");
	f_distribution.precision(N_PRECISION);
	f_distribution << "i " << "r_i" << endl;
	for (int i = 0; i < NUM_ITERATIONS; i++) {
		f_distribution << i << ' ' << drand48() << endl;
	}
	f_distribution.close();
	// Evaluation of distribution moments.
	f_moments.open("moments.dat");
	f_moments.precision(N_PRECISION);
	f_moments << "N " << "k=1 " << "k=2" << endl;
	N = 1;
	while (N <= MAX_ITERATIONS) {
		k_1 = 0.0;
		k_2 = 0.0;
		for (int i = 1; i <= N; i++) {
			x_i = drand48();
			k_1 += x_i;
			k_2 += x_i*x_i;
		}
		k_1 = fabs(k_1 / N - 0.5);
		k_2 = fabs(k_2 / N - 1.0 / 3.0);
		f_moments << N << ' ' << k_1 << ' ' << k_2 << endl;
		N *= 2;
	}
	f_moments.close();
	return 0;
}
