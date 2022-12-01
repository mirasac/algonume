#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

int main() {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);

	cout << "Practice session #1" << endl;
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

	cout << "\nPractice session #2" << endl;
	int const M = 4;
	double c[M] = {5, 16, 22, 15};
	A = mat_new(M, M);
	A[0][0] = 1;
	A[0][1] = 2;
	A[0][2] = 1;
	A[0][3] = -1;
	A[1][0] = 3;
	A[1][1] = 6;
	A[1][2] = 4;
	A[1][3] = 4;
	A[2][0] = 4;
	A[2][1] = 4;
	A[2][2] = 3;
	A[2][3] = 4;
	A[3][0] = 2;
	A[3][1] = 0;
	A[3][2] = 1;
	A[3][3] = 5;
	x = mat_new(1, M);
	x[0] = gaussian_elimination(A, c, M);
	cout << "Expected x = (4, -12, 22, -3)" << endl;
	cout << "Obtained x =" << endl;
	mat_cout(x, 1, M);
	mat_delete(A);
	mat_delete(x);
	return 0;
}

