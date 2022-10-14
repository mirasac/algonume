/*
In tis exercise we use the same number of gaussian points and the same interval for each dimension.
*/

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

double function2d(double x, double y) {
	return x*x*x*x * y*y + 2.0 * y*y * x*x - y * x*x + 2.0;
}

double function2d2(double x, double y) {
	double _return, g;
	g = sqrt(x*x + y*y);
	if (g <= 1.0) {
		_return = 1.0;
	} else {
		_return = 0.0;
	}
	return _return;
}

double simpsonquad(double (*f)(double x), double a, double b, int N);
double gaussquad(double (*f)(double x), double a, double b, int N, int Ng);
double multiquad(double (*f)(double x, double y), double a, double b, int N, int Ng);

int testquad(double (*q)(double (*f)(double x), double a, double b, int N), double (*f)(double x), double a, double b, double tollerance);

int main() {
	// Setup.
	using namespace std;
	int N, Ng;
	double a, b, tollerance, I, err, tmp;
	a = -1.0;
	b = 1.0;
	
	// Evaluate integral 1.
	cout << "Expected value: " << 412.0 / 45.0 << endl;
	cout << "Integral: " << multiquad(function2d, a, b, 1, 3) << endl;
	
	// Evaluate integral 2.
	// This example shows how gaussian quadrature is not good for discontinuous functions.
	cout << "Expected value: " << M_PI << endl;
	cout << "Integral: " << multiquad(function2d2, a, b, 1, 4) << endl;
	
	// Second part.
	N = 1;
	I = 0.0;
	tollerance = 1e-5;
	do {
		++N;
		tmp = multiquad(function2d2, a, b, N, 4);
		err = fabs(I - tmp);
		I = tmp;
	} while (err <= tollerance);
	cout << "Are needed N = " << N << " intervals" << endl;
	
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

int testquad(double (*q)(double (*f)(double x), double a, double b, int N), double (*f)(double x), double a, double b, double tollerance) {
	orderinterval(&a, &b);
	int N = 2;
	double t_0, t, err;
	t_0 = q(f, a, b, N);
	do {
		N *= 2;
		t = q(f, a, b, N);
		err = fabs(t - t_0);
		t_0 = t;
	} while (err > tollerance);
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

double multiquad(double (*f)(double x, double y), double a, double b, int N, int Ng) {
	orderinterval(&a, &b);
	double dx, dy, x_i, y_i, a_n, b_n, r_n, sum, sum_x, sum_y, sum_n;
	double * t, * w;
	int n_x, n_y;
	t = new double[Ng];
	w = new double[Ng];
	gausspoints(Ng, t, w);
	dx = (b - a) / N;
	dy = (b - a) / N;
	sum = 0.0;
	for (int i_x = 0; i_x < N; i_x++) {
		for (int i_y = 0; i_y < N; i_y++) {
			x_i = a + (i_x + 0.5) * dx;
			y_i = a + (i_y + 0.5) * dy;
			sum_y = 0.0;
			for (n_y = 0; n_y < Ng; n_y++) {
				sum_x = 0.0;
				for (n_x = 0; n_x < Ng; n_x++) {
					sum_x += w[n_x] * f(x_i + dx * t[n_x] / 2.0, y_i + dy * t[n_y] / 2.0);
				}
				sum_y += dx * w[n_y] * sum_x / 2.0;
			}
			sum += dy * sum_y / 2.0;
		}
	}
	delete[] t;
	delete[] w;
	return sum;
}

