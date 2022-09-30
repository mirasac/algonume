#include <iostream>

float add_one(float);

int main() {
	using namespace std;
	float a;
	cout << "Insert a number: ";
	cin >> a;
	cout << "Result: " << add_one(a) << endl;
	return 0;
}

float add_one(float a) {
	return a + 1.0;
}

