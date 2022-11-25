#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

int main() {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	int const N = 5;
	double ** A;
	double b[N] = {1, 0, 3, 1, 0};
	A = mat_new(N, N);
	mat_zero(A, N, N);
	A[0][0] = 2;
	A[1][1] = 2;
	A[2][2] = 2;
	A[3][3] = 2;
	A[4][4] = 2;
	A[1][0] = 1;
	A[2][1] = 1;
	A[3][2] = 1;
	A[4][3] = 1;
	A[0][1] = 1;
	A[1][2] = 1;
	A[2][3] = 1;
	A[3][4] = 1;
	cout << "Expected x = (2, -3, 4, -2, 1)" << endl;
	cout << "Result x =" << endl;
	// MC print result.
	return 0;
}

