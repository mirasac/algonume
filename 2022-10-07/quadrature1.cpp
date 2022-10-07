#include <iostream>
#include <cmath>
#include <cstdlib>

/*
Useful function I introduced to order the interval using the notation a inferior and b superior extremes.
*/
void orderinterval(double * a, double * b);

/*
Mathematical function I can change without editing the code for the quadrature.
*/
double function(double x) {
	return exp(-x);
}

double rectangualquad(double (*f)(double x), double a, double b, int N);
double midpointquad(double (*f)(double x), double a, double b, int N);
double trapezioidalquad(double (*f)(double x), double a, double b, int N);
double simpsonquad(double (*f)(double x), double a, double b, int N);

// MC testare.
int testquad(double (*q)(double (*f)(double x), double a, double b, int N), double (*f)(double x), double a, double b, double tollerance);

int main() {
	// Setup.
	using namespace std;
	int N;
	double a, b, tollerance, I;
	a = 0.0;
	b = 1.0;
	tollerance = 1e-9;
	
	// Test rectangular quadrature.
	N = testquad(rectangualquad, function, a, b, tollerance);
	I = rectangualquad(function, a, b, N);
	cout << "Rectangular: " << I << "; N = " << N << endl;
	
	// Test midpoint quadrature.
	N = testquad(midpointquad, function, a, b, tollerance);
	I = midpointquad(function, a, b, N);
	cout << "Midpoint: " << I << "; N = " << N << endl;
	
	// Test trapezioidal quadrature.
	N = testquad(trapezioidalquad, function, a, b, tollerance);
	I = trapezioidalquad(function, a, b, N);
	cout << "Trapezioidal: " << I << "; N = " << N << endl;
	
	// Test simpson quadrature.
	N = testquad(simpsonquad, function, a, b, tollerance);
	I = simpsonquad(function, a, b, N);
	cout << "Simpson: " << I << "; N = " << N << endl;
	return 0;
}

void orderinterval(double * a, double * b) {
	double * tmp;
	if (*a >= *b) {
		tmp = b;
		b = tmp;
		a = tmp;
	}
}

double rectangualquad(double (*f)(double x), double a, double b, int N) {
	orderinterval(&a, &b);
	double dx, x_n, s_n;
	dx = (b - a) / N;
	s_n = 0.0;
	// My implementation.
	/*
	x_n = a;
	for (int n = 1; n <= N: n++) {
		s_n += f(x_n) * dx;
		x_n += dx;
	}
	return s_n;
	*/
	// Better implementation because there are less repeated operations involved.
	for (int n = 0; n < N; n++) {
		x_n = a + n * dx;
		s_n += f(x_n);
	}
	return s_n * dx;
}

double midpointquad(double (*f)(double x), double a, double b, int N) {
	orderinterval(&a, &b);
	double dx, x_n, s_n;
	dx = (b - a) / N;
	x_n = a + dx / 2.0;
	s_n = 0.0;
	for (int n = 1; n <= N; n++) {
		s_n += f(x_n) * dx;
		x_n += dx;
	}
	return s_n;
}

double trapezioidalquad(double (*f)(double x), double a, double b, int N) {
	orderinterval(&a, &b);
	double dx, x_n, s_n;
	dx = (b - a) / N;
	s_n = 0.0;
	// My implementation, not wrong but inefficient because function f is called twice.
	x_n = a;
	for (int n = 1; n <= N; n++) {
		s_n += (f(x_n) + f(x_n + dx)) / 2.0 * dx;
		x_n += dx;
	}
	return s_n;
	// Better implementation, f is called only N + 1 times.
	// MC finire.
	/*
	for (int n = 0; n < N; n++) {
		x_n = a + n * dx;
		s_n = 
	}
	return s_n * dx;
	*/
}

// MC finire.
double simpsonquad(double (*f)(double x), double a, double b, int N) {
	orderinterval(&a, &b);
	if (N % 2 != 0) {
		std::cout << "Error: for the Simpson rule the number of intervals must be even" << std::endl;
		exit(1);
	}
	double dx, x_n, s_n;
	dx = (b - a) / N;
	x_n = a;
	s_n = (f(a) + f(b)) / 3.0;
	for (int n = 1; n < N; n++) {
		x_n += n * dx;
		s_n = (2.0 * f(x_n) + f(x_n + dx)) * 2.0 / 3.0;
	}
	return s_n * dx;
}

int testquad(double (*q)(double (*f)(double x), double a, double b, int N), double (*f)(double x), double a, double b, double tollerance) {
	orderinterval(&a, &b);
	int N = 1;
	double err;
	do {
		err = fabs(q(f, a, b, 2*N) - q(f, a, b, N));
		N *= 2;
	} while (err > tollerance);
	return N;
}

