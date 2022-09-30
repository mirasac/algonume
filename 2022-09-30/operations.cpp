#include <iostream>

int main() {
	using namespace std;
	float f1, f2;
	int i1, i2;
	cout << "Insert first float number: ";
	cin >> f1;
	cout << "Insert second float number: ";
	cin >> f2;
	cout << "Insert first int number: ";
	cin >> i1;
	cout << "Insert second int number: ";
	cin >> i2;
	cout << "Sum between floats: " << f1 + f2 << endl;
	cout << "Sum between ints: " << i1 + i2 << endl;
	cout << "Difference between floats: " << f1 - f2 << endl;
	cout << "Difference between ints: " << i1 - i2 << endl;
	cout << "Product between floats: " << f1 * f2 << endl;
	cout << "Product between ints: " << i1 * i2 << endl;
	cout << "Division between floats: " << f1 / f2 << endl;
	cout << "Division between ints: " << i1 / i2 << endl;
}

