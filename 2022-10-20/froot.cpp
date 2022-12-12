#include <cmath>
#include "../mclib/mclib.h"

#define ORDER 4
#define NUM_INTERVALS_5 10

double function(double x) {
	return exp(-x) - x;
}

double function_1(double x) {
	return -exp(-x) - 1.0;
}

double function2(double x) {
	return exp(1.0 / (x + 0.5)) - (3.0 + 2.0 * x) / (1.0 + x);
}

double function2_1(double x) {
	return - x / ( (x + 0.5) * (x + 0.5) ) * exp(1.0 / (x + 0.5)) + 1 / (1 + x*x);
}

double function3(double x) {
	return sin(x) - ( (x / 10.0) * (x / 10.0) + x / 5.0 + 1.0 / 3.0 );
}

int main() {
	using namespace std;
	double a, b, tolerance;

	cout << "Practice session #1" << endl;
	a = -1.0;
	b = 1.0;
	tolerance = 1e-7;
	bisection(function, a, b, tolerance);
	falseposition(function, a, b, tolerance);

	cout << endl << "Practice session #2" << endl;
	secant(function, a, b, tolerance);
	newtonraphson(function, function_1, a, b, tolerance);
	
	cout << endl << "Practice session #3" << endl;
	double c[ORDER];
	a = -5.0;
	b = 0.0;
	tolerance = 1e-8;
	c[3] = 1.0;
	c[2] = -3.0;
	c[1] = 1.0;
	c[0] = 5.0;
	newtonraphson_poly(polynomial, ORDER, c, a, b, tolerance);
	
	cout << endl << "Practice session #4" << endl;
	a = 0.0;
	b = 2.0;
	tolerance = 1e-7;
	bisection(function2, a, b, tolerance);
	falseposition(function2, a, b, tolerance);
	secant(function2, a, b, tolerance);
	newtonraphson(function2, function2_1, a, b, tolerance);
	
	cout << endl << "Practice session #5" << endl;
	a = -10.0;
	b = 10.0;
	double dx;
	double x_L[NUM_INTERVALS_5], x_R[NUM_INTERVALS_5];
	int N;
	N = bracketing(function3, a, b, NUM_INTERVALS_5, x_L, x_R);
	for (int n = 0; n < N; n++) {
		cout << "Root #" << n + 1 << "(bisection): " << bisection(function3, x_L[n], x_R[n], tolerance) << endl;
	}
}

