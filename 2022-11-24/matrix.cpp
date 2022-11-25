#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

int main() {
	using namespace std;
	cout << fixed << setprecision(1);
	int const N = 4;
	
	// Native method.
	cout << "Native method" << endl;
	double m_nat[N][N] = { {11, 12, 13, 14}, {21, 22, 23, 24}, {31, 32, 33, 34}, {41, 42, 43, 44} };
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << m_nat[i][j] << ' ';
		}
		cout << endl;
	}
	
	// Dynamical method.
	cout << "\nDynamical method" << endl;
	double ** m_dyn;
	m_dyn = mat_new(N, N);
	mat_zero(m_dyn, N, N);
	mat_cout(m_dyn, N, N);
	mat_delete(m_dyn);

	// Matrix-vector multiplication.	
	cout << "\nMatrix-vector multiplication" << endl;
	double ** A, ** b, ** Ab;
	A = mat_new(N, N);
	b = mat_new(N, 1);
	Ab = mat_new(N, 1);
	A[0][0] = 1;
	A[0][1] = 3;
	A[0][2] = 2;
	A[0][3] = -4;
	A[1][0] = 7;
	A[1][1] = 2;
	A[1][2] = 4;
	A[1][3] = 1;
	A[2][0] = 0;
	A[2][1] = -1;
	A[2][2] = 2;
	A[2][3] = 2;
	A[3][0] = 6;
	A[3][1] = 3;
	A[3][2] = 0;
	A[3][3] = 1;
	b[0][0] = 1;
	b[0][1] = 0;
	b[0][2] = 3;
	b[0][3] = 2;
	Ab = mat_multiply(A, b, N, N, 1);
	mat_cout(Ab, N, 1);
	mat_delete(A);
	mat_delete(b);
	mat_delete(Ab);
	return 0;
}

