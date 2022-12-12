#include <iostream>
#include <cmath>
#include <cstdlib>

/*
Useful function I introduced to order the interval using the notation a inferior and b superior extremes.
*/
void orderinterval(double * a, double * b);
void gausspoints(int Ng, double * x, double * w);

/*
Mathematical function I can change without editing the code for the quadrature.
*/
double function(double x) {
	return sqrt(1 + x);
}

double function2(double x) {
	return 1.0 - x + 2.0 * x * x + x * x * x / 2.0 + x * x * x * x / 4.0 - x * x * x * x * x / 8.0;
}

double simpsonquad(double (*f)(double x), double a, double b, int N);
double gaussquad(double (*f)(double x), double a, double b, int N, int Ng);

int testquad(double (*q)(double (*f)(double x), double a, double b, int N), double (*f)(double x), double a, double b, double tolerance);

int main() {
	// Setup.
	using namespace std;
	int N, Ng;
	double a, b, tolerance, I;
	a = 0.0;
	b = 3.0;
	
	// Evaluate integral 1.
	cout << "Integral 1" << endl;
	I = simpsonquad(function, a, b, 2);
	cout << "Simpson = " << I << endl;
	I = gaussquad(function, a, b, 1, 3);
	cout << "Gauss = " << I << endl;
	
	// Evaluat integral 2.
	// This example shows that gaussian quadrature can eavluat exactly the integral of a polynomial up to degree 2Ng - 1.
	cout << "Integral 2" << endl;
	a = -1.0;
	b = 5.0;
	cout << "Expected value: " << -66.0 / 5.0 << endl;
	I = simpsonquad(function2, a, b, 2);
	cout << "Simpson = " << I << endl;
	I = gaussquad(function2, a, b, 1, 3);
	cout << "Gauss = " << I << endl;
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

void gausspoints(int Ng, double * x, double * w) {
	switch (Ng) {
	case 1:
		x[0] = 0.0;
		w[0] = 2.0;
		break;
	case 2:
		x[0] = - sqrt(1.0 / 3.0);
		w[0] = 1.0;
		x[1] = sqrt(1.0 / 3.0);
		w[1] = 1.0;
		break;
	case 3:
		x[0] = - sqrt(3.0 / 5.0);
		w[0] = 5.0 / 9.0;
		x[1] = 0.0;
		w[1] = 8.0 / 9.0;
		x[2] = sqrt(3.0 / 5.0);
		w[2] = 5.0 / 9.0;
		break;
	// MC casi non utile per esercizio.
	case 4:
		std::cout << "Error: not yet implemented";
		exit(2);
		break;
	case 5:
		std::cout << "Error: not yet implemented";
		exit(2);
		break;
	}
}

double simpsonquad(double (*f)(double x), double a, double b, int N) {
	orderinterval(&a, &b);
	if (N % 2 != 0) {
		std::cout << "Error: to apply the Simpson rule the number of intervals must be even" << std::endl;
		exit(1);
	}
	double dx, x_n, s_n;
	dx = (b - a) / N;
	//x_n = a + dx;
	s_n = f(a) + f(b);
	for (int n = 1; n <= N - 1; n++) {
		x_n = a + n * dx;
		s_n += f(x_n) * 2.0 * (1 + n % 2);
	}
	return s_n * dx / 3.0;
}

int testquad(double (*q)(double (*f)(double x), double a, double b, int N), double (*f)(double x), double a, double b, double tolerance) {
	orderinterval(&a, &b);
	int N = 2;
	double t_0, t, err;
	t_0 = q(f, a, b, N);
	do {
		N *= 2;
		t = q(f, a, b, N);
		err = fabs(t - t_0);
		t_0 = t;
	} while (err > tolerance);
	return N;
}

double gaussquad(double (*f)(double x), double a, double b, int N, int Ng) {
	orderinterval(&a, &b);
	// MC usare malloc o new per creare array con Ng parametro.
	double dx, a_n, b_n, r_n, sum, sum_n;
	double * x, * w;
	x = new double[Ng];
	w = new double[Ng];
	gausspoints(Ng, x, w);
	dx = (b - a) / N;
	sum = 0.0;
	for (int n = 1; n <= N; n++) {
		sum_n = 0.0;
		a_n = a + (n - 1) * dx;
		b_n = a + n * dx;
		r_n = (b_n - a_n) / 2.0;
		for (int i = 0; i < Ng; i++) {
			sum_n += w[i] * f(r_n * x[i] + (b_n + a_n) / 2.0); 
		}
		sum_n = r_n * sum_n;
		sum += sum_n;
	}
	// MC capire bene come funziona new - delete.
	delete[] x;
	delete[] w;
	return sum;
}

