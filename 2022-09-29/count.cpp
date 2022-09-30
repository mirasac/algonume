#include <iostream>

#define MAX 10

int main() {
	using namespace std;
	cout << "Print with for loop:";
	int i; // Index declaration.
	for (i = 1; i <= MAX; i++) {
		cout << ' ' << i;
	}
	cout << endl;
	cout << "Print with while loop:";
	i = 1;
	while (i <= MAX) {
		cout << ' ' << i;
		i++;
	}
	cout << endl << "Print only odd numbers" << endl;
	cout << "Print with while loop:";
	for (i = 1; i <= MAX; i++) {
		if (i & 1) {
			cout << ' ' << i;
		}
	}
	cout << endl;
	cout << "Print with while loop:";
	i = 1;
	while (i <= MAX) {
		if (i & 1) {
			cout << ' ' << i;
		}
		i++;
	}
	cout << endl;
}

