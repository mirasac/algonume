#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

#define NX 5
#define NY NX
#define IBEG 1
#define IEND (NX - 2)
#define JBEG 1
#define JEND (NY - 2)

static double const global_S = 0.0;

double solution(double x, double y) {
	return exp(-M_PI * x) * sin(-M_PI * y) + global_S / 4.0 * (x*x + y*y);
}

void set_boundary_conditions(double ** phi, double x[], double y[]) {
	for (int i = 0; i < NX; i++) {
		phi[i][0] = solution(x[i], y[0]);
		phi[i][NY-1] = solution(x[i], y[NY-1]);
	}
	for (int j = 1; j < NY; j++) {
		phi[0][j] = solution(x[0], y[j]);
		phi[NX-1][j] = solution(x[NX-1], y[j]);
	}
}

double residual(double ** phi, double ** S, double h) {
	double epsilon, delta_x, delta_y;
	epsilon = 0.0;
	for (int i = 1; i <= NX - 2; i++) {
	for (int j = 1; j <= NY - 2; j++) {
		delta_x = phi[i+1][j] - 2.0 * phi[i][j] + phi[i-1][j];
		delta_y = phi[i][j+1] - 2.0 * phi[i][j] + phi[i][j-1];
		epsilon += abs(delta_x + delta_y - h*h * S[i][j]);
	}
	}
	return epsilon;
}

int main() {
	using namespace std;
	cout << setprecision(N_PRECISION);
	double const tolerance = 1e-7;
	double const x_min = 0.0;
	double const x_max = 1.0;
	double const y_min = 0.0;
	double const y_max = 1.0;
	int c;
	double dx, dy, epsilon;
	double ** phi_0, ** phi_1, ** S, x[NX], y[NY];
	dx = (x_max - x_min) / (NX - 1);
	dy = (y_max - y_min) / (NY - 1);
	phi_0 = mat_new(NX, NY);
	phi_1 = mat_new(NX, NY);
	S = mat_new(NX, NY);
	// Initialize internal points and assign boundary conditions.
	for (int i = 0; i < NX; i++) {
		x[i] = x_min + i * dx;
	}
	for (int j = 0; j < NY; j++) {
		y[j] = y_min + j * dy;
	}
	mat_zero(phi_0, NX, NY);
	mat_zero(phi_1, NX, NY);
	set_boundary_conditions(phi_0, x, y);
	mat_constant(S, global_S, NX, NY);
	//mat_cout(phi_0, NX, NY); // MC debug.
	// Find solution with Gauss-Seidel.
	c = 0;
	do {
		++c;
		gauss_seidel(phi_0, S, dx, NX, NY);
		mat_cout(phi_0, NX, NY); // MC debug.
		epsilon = residual(phi_0, S, dx);
	} while (epsilon > tolerance);
	cout << "Gauss-Siedel" << endl;
	cout << "S = " << global_S;
	cout << ", " << c << " iterations" << endl;
	// MC convergence after only 1 step: something does not work.
	// Clean up.
	mat_delete(phi_0);
	mat_delete(phi_1);
	mat_delete(S);
	return 0;
}
