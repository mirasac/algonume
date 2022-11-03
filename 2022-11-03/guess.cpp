#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

int main() {
	using namespace std;
	srand(time(NULL));
	cout << setprecision(N_PRECISION);
	int target, guess, a, b, c;
	a = 1;
	b = 100;
	c = 0;
	target = rand() % (b - a + 1) + a;
	while (guess != target) {
		cout << "n in [" << a << "," << b << "]" << endl;
		cout << "type your guess #" << ++c << endl;
		cin >> guess;
		if (guess == target) {
			cout << "Number found in " << c << " iterations !!" << endl;
		} else if (guess < target) {
			a = guess;
		} else if (guess > target) {
			b = guess;
		}
	}
	return 0;
}
