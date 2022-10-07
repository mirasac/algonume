/*
Use fabs for float absolute value.
*/

#include <cmath>
#include <iostream>
#include <iomanip>

#define PRECISION 12
#define TOLLERANCE 1e-9

float iteration(float x, float S) {
	return 0.5 * (x + S / x);
}

int main() {
	// Setup.
	using namespace std;
	cout << setiosflags(ios::scientific) << setprecision(PRECISION);
	float S, x, y, err;
	
	// Input data.
	cout << "Enter a real number: ";
	cin >> S;
	cout << "Enter your guess: ";
	cin >> x;
	cout << "---------" << endl;
	
	// Iterations.
	int n = 0;
	do {
		++n;
		cout << "Iteration #" << n << ": ";
		y = iteration (x, S);
		err = fabs(y - x);
		cout << "x = " << y << " err = " << err << endl;
		x = y;
	} while (err > TOLLERANCE);
	
	// Final messages.
	cout << "The SQRT of " << S << " is : " << x << endl;
	cout << "The exact values is: " << sqrt(S) << endl;
}

