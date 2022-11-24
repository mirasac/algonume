#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	int const N = 4;
	
	// Native method.
	double m[N][N];
	mat_zero(m, N);
	mat_print(m, N);
	
	// Dynamical method.
	//mat_print(m, N);
	
	return 0;
}

