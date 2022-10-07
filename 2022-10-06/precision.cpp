#include <iostream>

int main() {
	using namespace std;
	float epsilon_f = 1.0;
	while ((float) 1.0 + epsilon_f != (float) 1.0) {
		epsilon_f /= (float) 10.0;
	}
	cout << "Machine epsilon for float type is: " << epsilon_f << endl;
	double epsilon_d = 1.0;
	while ((double) 1.0 + epsilon_d != (double) 1.0) {
		epsilon_d /= (double) 10.0;
	}
	cout << "Machine epsilon for double type is: " << epsilon_d << endl;
	return 0;
}
