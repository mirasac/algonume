#include <iostream>
#include <cmath>

double eval_b(double x1, double x2) {
	return - (x1 + x2);
}

double eval_c(double x1, double x2) {
	return x1 * x2;
}

void solve_standard(double a, double b, double c, double &x_plus, double &x_minus) {
	x_plus = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
	x_minus = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
}

void solve_sign(double a, double b, double c, double &x_plus, double &x_minus) {
	if (b >= 0.0) {
		x_plus = 2 * c / (-b - sqrt(b * b - 4 * a * c));
		x_minus = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
	} else {
		x_plus = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
		x_minus = (2 * c) / (-b + sqrt(b * b - 4 * a * c));
	}
}

void test(double a, double b, double c) {
	using namespace std;
	double x_plus, x_minus;
	solve_standard(a, b, c, x_plus, x_minus);
	cout << "With standard solving" << endl;
	cout << "x_plus = " << x_plus << endl;
	cout << "x_minus = " << x_minus << endl;
	solve_sign(a, b, c, x_plus, x_minus);
	cout << "With solving considering sign of beta" << endl;
	cout << "x_plus = " << x_plus << endl;
	cout << "x_minus = " << x_minus << endl;
}

int main() {
	using namespace std;
	double a, b, c, x1, x2;

	// Test 1
	x1 = 2.0;
	x2 = -3.0;
	a = 1.0;
	b = eval_b(x1, x2);
	c = eval_c(x2, x2);
	test(a, b, c);
	
	// Test 2
	x1 = 1e-5;
	x2 = 1e8;
	a = 1.0;
	b = eval_b(x1, x2);
	c = eval_c(x1, x2);
	test(a, b, c);

	// Test 3
	x1 = 1e-12;
	x2 = 1e12;
	a = 1.0;
	b = eval_b(x1, x2);
	c = eval_c(x1, x2);
	test(a, b, c);

	return 0;
}
