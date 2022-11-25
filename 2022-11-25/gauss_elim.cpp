#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

int main() {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	int const N = 3;
	double ** A, ** x;
	double b[N] = {1, 2, 1};
	A = mat_new(N, N);
	A[0][0] = 2;
	A[0][1] = -1;
	A[0][2] = 0;
	A[1][0] = -1;
	A[1][1] = 2;
	A[1][2] = -1;
	A[2][0] = 0;
	A[2][1] = -1;
	A[2][2] = 2;
	x = mat_new(1, N);
	x[0] = gaussian_elimination(A, b, N);
	cout << "Expected x = (2, 3, 2)" << endl;
	cout << "Obtained x =" << endl;
	mat_cout(x, 1, N);
	mat_delete(A);
	mat_delete(x);
	return 0;
}

