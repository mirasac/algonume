#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

#define TIME_THRESHOLD 500

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	srand48(time(NULL));
	double lambda;
	int n_tot, n_dec, delta_t;
	ofstream f_plot;
	lambda = 0.01;
	n_tot = 1200000;
	delta_t = 1;
	f_plot.open("decay.dat");
	f_plot.precision(N_PRECISION);
	f_plot << "t " << "n_dec" << endl;
	for (int t = 1; t <= TIME_THRESHOLD || n_tot <= 0; t += delta_t) {
		f_plot << t << ' ';
		n_dec = 0;
		for (int n = 1; n <= n_tot; n++) {
			if (drand48() <= lambda) {
				++n_dec;
			}
		}
		f_plot << n_dec << endl;
		n_tot -= n_dec;
	}
	f_plot.close();
	return 0;
}

