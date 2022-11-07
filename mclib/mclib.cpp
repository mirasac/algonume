#include "mclib.h"



////////// Utilities functions //////////



void orderinterval(double * a, double * b) {
	double tmp;
	if (*a >= *b) {
		tmp = *b;
		*b = *a;
		*a = tmp;
	}
}

void gaussianpoints(int Ng, double t[], double w[]) {
	switch (Ng) {
	case 1:
		t[0] = 0.0;
		w[0] = 2.0;
		break;
	case 2:
		t[0] = -sqrt(1.0 / 3.0);
		w[0] = 1.0;
		t[1] = sqrt(1.0 / 3.0);
		w[1] = 1.0;
		break;
	case 3:
		t[0] = -sqrt(3.0 / 5.0);
		w[0] = 5.0 / 9.0;
		t[1] = 0.0;
		w[1] = 8.0 / 9.0;
		t[2] = sqrt(3.0 / 5.0);
		w[2] = 5.0 / 9.0;
		break;
	case 4:
		t[0] = -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
		w[0] = (18.0 + sqrt(30.0)) / 36.0;
		t[1] = -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
		w[1] = (18.0 - sqrt(30.0)) / 36.0;
		t[2] = sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
		w[2] = (18.0 + sqrt(30.0)) / 36.0;
		t[3] = sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
		w[3] = (18.0 - sqrt(30.0)) / 36.0;
		break;
	case 5:
		t[0] = -sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0;
		w[0] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
		t[1] = -sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0;
		w[1] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
		t[2] = 0.0;
		w[2] = 128.0 / 225.0;
		t[3] = sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0;
		w[3] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
		t[4] = sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0;
		w[4] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
		break;
	default:
		std::cout << "Error: coefficients for Ng = " << Ng << " not available";
		exit(1);
		break;
	}
}

char check_intermediate_value(double (*f)(double x), double a, double b) {
	char _return = 0;
	if (f(a) * f(b) <= 0.0) {
		_return = 1;
	}
	return _return;
}

double polynomial(int n, double c[], double x) {
	double _return;
	_return = c[n];
	for (int i = n-1; i >= 0; i--) {
		_return = _return * x + c[i];
	}
	return _return;
}

void polynomial_derivative(int n, double c[], double c1[]) {
	for (int i = 1; i <= n; i++) {
		c1[i-1] = c[i] * i;
	}
}

double polynomial_legendre(int n, double x) {
	double _return, P_h, P_i;
	P_h = 1.0; // Order 0.
	P_i = x; // Order 1.
	switch (n) {
	case 0:
		_return = P_h;
		break;
	case 1:
		_return = P_i;
		break;
	default:
		for (int i = 1; i < n; i++) {
			_return = ( (2 * i + 1) * x * P_i - i * P_h ) / (i + 1);
			P_h = P_i;
			P_i = _return;
		}
		break;
	}
	return _return;
}



////////// Quadrature functions //////////



double rectangularquad(double (*f)(double x), double a, double b, int N) {
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

double gaussquad(double (*f)(double x), double a, double b, int N, int Ng) {
	orderinterval(&a, &b);
	// MC usare malloc o new per creare array con Ng parametro.
	double dx, a_n, b_n, r_n, sum, sum_n;
	double * x, * w;
	x = new double[Ng];
	w = new double[Ng];
	gaussianpoints(Ng, x, w);
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
	gaussianpoints(Ng, t, w);
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



////////// Root finders //////////



// MC one can improve this function avoiding multiple calls to f by saving the values in a variable.
double bisection(double (*f)(double x), double a, double b, double tollerance) {
	if (!check_intermediate_value(f, a, b)) {
		std::cout << "Function f, assumed continuous, does not have zeros in interval [" << a << ", " << b << "]";
		exit(1);
	}
	double x_0;
	#if FLG_DEBUG
	int k = 0;
	#endif /* FLG_DEBUG */
	do {
		x_0 = (a + b) / 2.0;
		#if FLG_DEBUG
		++k;
		std::cout << "Bisection(): k = " << k << "; [a,b] = [" << a << ", " << b << "]; xm = " << x_0 << "; dx = " << b - a << "; fm = " << f(x_0) << std::endl;
		#endif /* FLG_DEBUG */
		if (check_intermediate_value(f, a, x_0)) {
			b = x_0;
		} else {
			a = x_0;
		}
	} while (fabs(b - a) > tollerance);
	#if FLG_DEBUG
	++k;
	std::cout << "Bisection(): k = " << k << "; [a,b] = [" << a << ", " << b << "]; xm = " << x_0 << "; dx = " << b - a << "; fm = " << f(x_0) << std::endl;
	#endif /* FLG_DEBUG */
	return x_0;
}

double falseposition(double (*f)(double x), double a, double b, double tollerance) {
	orderinterval(&a, &b);
	if (!check_intermediate_value(f, a, b)) {
		std::cout << "Function f, assumed continuous, does not have zeros in interval [" << a << ", " << b << "]";
		exit(1);
	}
	double x_m, f_a, f_b, f_m;
	#if FLG_DEBUG
	int k = 0;
	#endif /* FLG_DEBUG */
	f_a = f(a);
	f_b = f(b);
	while (b - a >= tollerance) {
		x_m = (f_a * b - f_b * a) / (f_a - f_b);
		f_m = f(x_m);
		#if FLG_DEBUG
		++k;
		std::cout << "FalsePosition(): k = " << k << "; [a,b] = [" << a << ", " << b << "]; xm = " << x_m << "; dx = " << b - a << "; fm = " << f_m << std::endl;
		#endif /* FLG_DEBUG */
		if (f_a * f_m >= 0.0) {
			a = x_m;
			f_a = f_m;
		} else {
			b = x_m;
			f_b = f_m;
		}
	}
	return x_m;
}

double secant(double (*f)(double x), double a, double b, double tollerance) {
	double x_0, x_prev, f_0, f_prev;
	#if FLG_DEBUG
	int k = 0;
	#endif /* FLG_DEBUG */
	x_0 = (a + b) / 2.0;  // Initial guess.
	x_prev = a;
	f_prev = f(x_prev);
	do {
		f_0 = f(x_0);
		x_0 = x_0 - f_0 * (x_0 - x_prev) / (f_0 - f_prev);
		x_prev = x_0;
		f_prev = f_0;
		#if FLG_DEBUG
		++k;
		std::cout << "Secant(): k = " << k << "; xc = " << x_0 << "; dx = " << x_0 - x_prev << std::endl;
		#endif /* FLG_DEBUG */
	} while(fabs(x_0 - x_prev) >= tollerance);
	#if FLG_DEBUG
	++k;
	std::cout << "Secant(): k = " << k << "; xc = " << x_0 << "; dx = " << x_prev - x_0 << std::endl;
	std::cout << "Secant(): f = " << f(x_0) << std::endl;
	#endif /* FLG_DEBUG */
	return x_0;
}

double newtonraphson(double (*f)(double x), double (*f1)(double x), double a, double b, double tollerance) {
	orderinterval(&a, &b);
	double x_prev, x_0, f_0;
	#if FLG_DEBUG
	int k = 0;
	#endif /* FLG_DEBUG */
	x_0 = (a + b) / 2.0;  // Initial guess.
	do {
		f_0 = f(x_0);
		x_prev = x_0;
		x_0 = x_0 - f_0 / f1(x_0);
		#if FLG_DEBUG
		++k;
		std::cout << "Newton(): k = " << k << "; xc = " << x_0 << "; dx = " << x_0 - x_prev << std::endl;
		#endif /* FLG_DEBUG */
	} while (fabs(x_0 - x_prev) >= tollerance);
	#if FLG_DEBUG
	++k;
	std::cout << "Newton(): k = " << k << "; xc = " << x_0 << "; dx = " << x_prev - x_0 << std::endl;
	std::cout << "Newton(): f = " << f(x_0) << std::endl;
	#endif /* FLG_DEBUG */
	return x_0;
}

double newtonraphson_poly(double (*p)(int n, double c[], double x), int n, double c[], double a, double b, double tollerance) {
	orderinterval(&a, &b);
	double x_prev, x_0, f_0, f1;
	double * c1;
	c1 = new double[n];
	polynomial_derivative(n, c, c1);
	#if FLG_DEBUG
	int k = 0;
	#endif /* FLG_DEBUG */
	x_0 = (a + b) / 2.0;  // Initial guess.
	do {
		f_0 = p(n, c, x_0);
		f1 = p(n-1, c1, x_0);
		x_prev = x_0;
		x_0 = x_0 - f_0 / f1;
		#if FLG_DEBUG
		++k;
		std::cout << "Newton(): k = " << k << "; xc = " << x_0 << "; dx = " << x_0 - x_prev << std::endl;
		#endif /* FLG_DEBUG */
	} while (fabs(x_0 - x_prev) >= tollerance);
	#if FLG_DEBUG
	++k;
	std::cout << "Newton(): k = " << k << "; xc = " << x_0 << "; dx = " << x_prev - x_0 << std::endl;
	std::cout << "Newton(): f = " << p(n, c, x_0) << std::endl;
	#endif /* FLG_DEBUG */
	delete[] c1;
	return x_0;
}



////////// Bracketing //////////



int bracketing(double (*f)(double x), double a, double b, int N, double x_L[], double x_R[]) {
	orderinterval(&a, &b);
	double x_a, x_b, dx;
	int k;
	dx = (b - a) / N;
	k = 0;
	for (int n = 0; n < N; n++) {
		x_a = a + dx * n;
		x_b = a + dx * (n + 1);
		if (check_intermediate_value(f, x_a, x_b)) {
			x_L[k] = x_a;
			x_R[k] = x_b;
			++k;
		}
	}
	return k;
}



////////// Differentation //////////



double forward_difference(double (*f)(double), double h, double x) {
	return (f(x + h) - f(x)) / h;
}

double backward_difference(double (*f)(double), double h, double x) {
	return (f(x) - f(x - h)) / h;
}

double central_difference(double (*f)(double), double h, double x) {
	return (f(x + h) - f(x - h)) / (2.0 * h);
}



////////// ODE integration //////////



// MC finire.
/*
void euler_step(void (*RHS)(double t, double * y_0, double * y), double t, double * y, double dt, int n_eq) {
	int k;
	double * rhs;
	rhs = new double[n_eq];
	RHS(t, 
}
*/
