#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

int main() {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	int const N = 5;
	double * x;
	double d_inf[N] = {0.0, 1.0, 1.0, 1.0, 1.0};
	double d[N] = {2.0, 2.0, 2.0, 2.0, 2.0};
	double d_sup[N] = {1.0, 1.0, 1.0, 1.0, 0.0};
	double b[N] = {1.0, 0.0, 3.0, 1.0, 0.0};
	/*double ** A;
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
	A[3][4] = 1;*/
	cout << "Expected x = (2, -3, 4, -2, 1)" << endl;
	cout << "Result x is:" << endl;
	x = tridiagonal_solver(d_inf, d, d_sup, b, N);
	for (int i = 0; i < N; i++) {
		cout << x[i] << endl;
	}
	delete[] x;
	return 0;
}

