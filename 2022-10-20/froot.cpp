#include <cmath>
#include "../mclib/mclib.h"

#define ORDER 4

double function(double x) {
	return exp(-x) - x;
}

double function1(double x) {
	return -exp(-x) - 1.0;
}

double function2(double x) {
	return exp(1.0 / (x + 0.5)) - (3.0 + 2.0 * x) / (1.0 + x);
}

int main() {
	using namespace std;
	double a, b, tollerance;

	cout << "Practice session #1" << endl;
	a = -1.0;
	b = 1.0;
	tollerance = 1e-7;
	bisection(function, a, b, tollerance);
	falseposition(function, a, b, tollerance);

	cout << endl << "Practice session #2" << endl;
	secant(function, a, b, tollerance);
	newtonraphson(function, function1, a, b, tollerance);
	
	cout << endl << "Practice session #3" << endl;
	double c[ORDER];
	a = -5.0;
	b = 0.0;
	tollerance = 1e-8;
	c[3] = 1.0;
	c[2] = -3.0;
	c[1] = 1.0;
	c[0] = 5.0;
	newtonraphson_poly(polynomial, ORDER, c, a, b, tollerance);
	
	cout << endl << "Practice session #4" << endl;
	a = 0.0;
	b = 2.0;
	tollerance = 1e-7;
	secant(function2, a, b, tollerance);
}

