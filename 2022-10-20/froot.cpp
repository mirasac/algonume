#include <cmath>
#include "../mclib/mclib.h"

double function(double x) {
	return exp(-x) - x;
}

int main() {
	double a, b, tollerance;
	a = -1.0;
	b = 1.0;
	tollerance = 1e-7;
	bisection(function, a, b, tollerance);
}
