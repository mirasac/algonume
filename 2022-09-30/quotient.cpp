#include <iostream>

int quotient(int, int, float &, float &);

int main() {
	using namespace std;
	int a, b;
	cout << "Insert numerator: ";
	cin >> a;
	cout << endl << "Insert denominator: ";
	cin >> b;
	cout << endl;
	float q, r;
	if (quotient(a, b, q, r)) {
		cout << "Quotient: " << q << endl;
		cout << "Remainder: " << r << endl;
	} else {
		cout << "Error: denominator is 0";
	}
	return 0;
}

int quotient(int n, int d, float &q, float &r) {
	if (d != 0) {
		q = n / d;
		r = n % d;
		return 1;
	} else {
		return 0;
	}
}

