#include <iostream>

int sum_int(int, int);

float sum_float(float, float);

int main() {
	using namespace std;
	float a, b;
	cout << "Insert first number: ";
	cin >> a;
	cout << "Insert second number: ";
	cin >> b;
	cout << "Sum int: " << sum_int(a, b) << endl;
	cout << "Sum float: " << sum_float(a, b) << endl;
	return 0;
}

int sum_int(int a, int b) {
	return a + b;
}

float sum_float(float a, float b) {
	return a + b;
}

