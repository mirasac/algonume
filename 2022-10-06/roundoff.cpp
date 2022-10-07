/*
Get numerical approx sqrt(x^2 + 1) = x x >> 1 and 1 - cos(x) x~0
*/

#include <iostream>
#include <cmath>

float fx1_sqrt(float x) {
	return sqrt(x * x + 1.0) - x;
}

float fx2_sqrt(float x) {
	return 1.0 / (sqrt(x * x + 1.0) + x);
}

float taylor_sqrt(float x) {
	return 1.0 / (sqrt(x * x + 1.0) + x); // TODO
}

float fx1_cos(float x) {
	return 1.0 - cos(x);
}

float fx2_cos(float x) {
	return 1.0 - cos(x); // How to rationalize cos? Is it Taylor again? Or is it discretized with steps large eps?
}

float taylor_cos(float x) {
	return 1.0 - cos(x); // TODO
}

int main() {
	using namespace std;
	
	//cout << setiosflags(ios::scientific); // Another library is needed.
	
	cout << "Example #1: compute sqrt(x^2 + 1) ~ x for large x" << endl;
	cout << "---------" << endl;
	// Second example.
	cout << "Example #2: compute sqrt(x^2 + 1) ~ x for large x" << endl;
	cout << "---------" << endl;
}

