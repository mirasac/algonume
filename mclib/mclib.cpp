#include "mclib.h"



////////// Utilities functions //////////



void orderinterval(double * a, double * b) {
	double * tmp;
	if (*a >= *b) {
		tmp = b;
		b = tmp;
		a = tmp;
	}
}

void gausspoints(int Ng, double * x, double * w) {
	switch (Ng) {
	case 1:
		x[0] = 0.0;
		w[0] = 2.0;
		break;
	case 2:
		x[0] = -sqrt(1.0 / 3.0);
		w[0] = 1.0;
		x[1] = sqrt(1.0 / 3.0);
		w[1] = 1.0;
		break;
	case 3:
		x[0] = -sqrt(3.0 / 5.0);
		w[0] = 5.0 / 9.0;
		x[1] = 0.0;
		w[1] = 8.0 / 9.0;
		x[2] = sqrt(3.0 / 5.0);
		w[2] = 5.0 / 9.0;
		break;
	case 4:
		x[0] = -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
		w[0] = (18.0 + sqrt(30.0)) / 36.0;
		x[1] = -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
		w[1] = (18.0 - sqrt(30.0)) / 36.0;
		x[2] = sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
		w[2] = (18.0 + sqrt(30.0)) / 36.0;
		x[3] = sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
		w[3] = (18.0 - sqrt(30.0)) / 36.0;
		break;
	case 5:
		x[0] = -sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0;
		w[0] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
		x[1] = -sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0;
		w[1] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
		x[2] = 0.0;
		w[2] = 128.0 / 225.0;
		x[3] = sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0;
		w[3] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
		x[4] = sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0;
		w[4] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
		break;
	default:
		std::cout << "Error: coefficients for Ng = " << Ng << " not available";
		exit(1);
		break;
	}
}


////////// Quadrature functions //////////



