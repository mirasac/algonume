#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

#define NX 128
#define NY NX
#define IBEG 1
#define IEND (NX - 2)
#define JBEG 1
#define JEND (NY - 2)
#define PRACTICE_SESSION 2

static double const global_S = 0.0;

double solution1(double x, double y) {
	return exp(-M_PI * x) * sin(-M_PI * y) + global_S / 4.0 * (x*x + y*y);
}

// Only Dirichlet boundary conditions can be assigned only the first time, Neumann boundary conditions must be applied at each iteration.
void set_boundary_conditions1(double ** phi, double x[], double y[]) {
	for (int i = 0; i < NX; i++) {
		phi[i][0] = solution1(x[i], y[0]);
		phi[i][NY-1] = solution1(x[i], y[NY-1]);
	}
	for (int j = 1; j < NY - 1; j++) {
		phi[0][j] = solution1(x[0], y[j]);
		phi[NX-1][j] = solution1(x[NX-1], y[j]);
	}
}

static double const global_rho_0 = 1.0;
static double const global_a = 0.1;

double rho(double x, double y) {
	double r, _return;
	r = sqrt(x*x + y*y);
	if (r <= a) {
		_return = global_rho_0;
	} else {
		_return = 0.0;
	}
	return _return;
}

double solution2(double x, double y) {
	double r, _return;
	r = sqrt(x*x + y*y);
	if (r >= 0 && r <= a) {
		_return = - global_rho_0 * r*r / 4.0;
	} else {
		_return = - global_rho_0 * global_a*global_a / 2.0 * (log(r / a) + 0.5);
	}
	return _return;
}

void set_boundary_conditions2(double ** phi, double x[], double y[]) {
	for (int i = 0; i < NX; i++) {
		phi[i][0] = solution2(x[i], y[0]);
		phi[i][NY-1] = solution2(x[i], y[NY-1]);
	}
	for (int j = 1; j < NY - 1; j++) {
		phi[0][j] = solution2(x[0], y[j]);
		phi[NX-1][j] = solution2(x[NX-1], y[j]);
	}
}

double residual(double ** phi, double ** S, double dx, double dy) {
	double epsilon, delta_x, delta_y;
	epsilon = 0.0;
	for (int i = IBEG; i <= IEND; i++) {
	for (int j = JBEG; j <= JEND; j++) {
		delta_x = phi[i+1][j] - 2.0 * phi[i][j] + phi[i-1][j];
		delta_y = phi[i][j+1] - 2.0 * phi[i][j] + phi[i][j-1];
		epsilon += fabs(delta_x + delta_y - dx*dy * S[i][j]);
	}
	}
	return epsilon;
}

int main() {
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	double const tolerance = 1e-7;
	int c;
	double x_min, x_max, y_min, y_max, dx, dy, epsilon, omega;
	double ** phi, ** S, x[NX], y[NY];
	ofstream plot_file;
	
	
	
	// Practice session 1.
	#if PRACTICE_SESSION == 1
	cout << "Practice session #1" << endl;
	cout << "S = " << global_S << endl;
	x_min = 0.0;
	x_max = 1.0;
	y_min = 0.0;
	y_max = 1.0;
	dx = (x_max - x_min) / (NX - 1);
	dy = (y_max - y_min) / (NY - 1);
	// Initialize grid.
	for (int i = 0; i < NX; i++) {
		x[i] = x_min + i * dx;
	}
	for (int j = 0; j < NY; j++) {
		y[j] = y_min + j * dy;
	}
	
	// Solve with Jacobi method.
	phi = mat_new(NX, NY);
	S = mat_new(NX, NY);
	// Assign internal points and assign boundary conditions.
	mat_zero(phi, NX, NY);
	set_boundary_conditions1(phi, x, y);
	mat_constant(S, global_S, NX, NY);
	// Find solution1.
	c = 0;
	epsilon = residual(phi, S, dx, dx);
	while (epsilon > tolerance) {
		++c;
		jacobi(phi, S, dx, NX, NY);
		epsilon = residual(phi, S, dx, dx);
	}
	--c;
	// Print results on file.
	plot_file.open("elliptic_1_jacobi.dat");
	plot_file << fixed << setprecision(N_PRECISION);
	plot_file << "x y phi(x,y)" << endl;
	for (int j = 0; j < NY; j++) {
		for (int i = 0; i < NX; i++) {
			plot_file << x[i] << ' ' << y[j] << ' ' << phi[i][j] << endl;
		}
		plot_file << endl;
	}
	cout << "Jacobi, solved in " << c << " iterations, file printed" << endl;
	// Clean up.
	plot_file.close();
	mat_delete(phi);
	mat_delete(S);
	
	// Solve with Gauss-Seidel method.
	phi = mat_new(NX, NY);
	S = mat_new(NX, NY);
	// Assign internal points and assign boundary conditions.
	mat_zero(phi, NX, NY);
	set_boundary_conditions1(phi, x, y);
	mat_constant(S, global_S, NX, NY);
	// Find solution1.
	c = 0;
	epsilon = residual(phi, S, dx, dx);
	while (epsilon > tolerance) {
		++c;
		gauss_seidel(phi, S, dx, NX, NY);
		epsilon = residual(phi, S, dx, dx);
	}
	--c;
	// Print results on file.
	plot_file.open("elliptic_1_gauss_seidel.dat");
	plot_file << fixed << setprecision(N_PRECISION);
	plot_file << "x y phi(x,y)" << endl;
	for (int j = 0; j < NY; j++) {
		for (int i = 0; i < NX; i++) {
			plot_file << x[i] << ' ' << y[j] << ' ' << phi[i][j] << endl;
		}
		plot_file << endl;
	}
	cout << "Gauss-Seidel, solved in " << c << " iterations, file printed" << endl;
	// Clean up.
	plot_file.close();
	mat_delete(phi);
	mat_delete(S);
	
	// Solve with successive over relaxation method.
	omega = 2.0 / (1.0 + M_PI / NX);
	phi = mat_new(NX, NY);
	S = mat_new(NX, NY);
	// Assign internal points and assign boundary conditions.
	mat_zero(phi, NX, NY);
	set_boundary_conditions1(phi, x, y);
	mat_constant(S, global_S, NX, NY);
	// Find solution.
	c = 0;
	epsilon = residual(phi, S, dx, dx);
	while (epsilon > tolerance) {
		++c;
		successive_over_relaxation(phi, S, dx, omega, NX, NY);
		epsilon = residual(phi, S, dx, dx);
	}
	--c;
	// Print results on file.
	plot_file.open("elliptic_1_SOR.dat");
	plot_file << fixed << setprecision(N_PRECISION);
	plot_file << "x y phi(x,y)" << endl;
	for (int j = 0; j < NY; j++) {
		for (int i = 0; i < NX; i++) {
			plot_file << x[i] << ' ' << y[j] << ' ' << phi[i][j] << endl;
		}
		plot_file << endl;
	}
	cout << "SOR, solved in " << c << " iterations, file printed" << endl;
	// Clean up.
	plot_file.close();
	mat_delete(phi);
	mat_delete(S);
	#endif /* PRACTICE_SESSION */
	
	
	
	// Practice session 2.
	#if PRACTICE_SESSION == 2
	cout << "\nPractice session #2" << endl;
	x_min = -1.0;
	x_max = 1.0;
	y_min = -1.0;
	y_max = 1.0;
	dx = (x_max - x_min) / (NX - 1);
	dy = (y_max - y_min) / (NY - 1);
	// Initialize grid.
	for (int i = 0; i < NX; i++) {
		x[i] = x_min + i * dx;
	}
	for (int j = 0; j < NY; j++) {
		y[j] = y_min + j * dy;
	}
	// Jacobi.
	phi = mat_new(NX, NY);
	S = mat_new(NX, NY);
	// Assign initial values and boundary conditions.
	mat_zero(phi, NX, NY);
	set_boundary_conditions2(phi, x, y);
	for (int i = 0; i < NX, i++) {
		for (int j = 0; j < NY, j++) {
			S[i][j] = -rho(x[i], y[i]);
		}
	}
	// Find solution.
	c = 0;
	epsilon = residual(phi, S, dx, dx);
	while (epsilon > tolerance) {
		++c;
		jacobi(phi, S, dx, NX, NY);
		epsilon = residual(phi, S, dx, dx);
	}
	--c;
	// Print results on file.
	plot_file.open("elliptic_2_jacobi.dat");
	plot_file << fixed << setprecision(N_PRECISION);
	plot_file << "x y phi(x,y)" << endl;
	for (int j = 0; j < NY; j++) {
		for (int i = 0; i < NX; i++) {
			plot_file << x[i] << ' ' << y[j] << ' ' << phi[i][j] << endl;
		}
		plot_file << endl;
	}
	cout << "Jacobi solved in " << c << " iterations, file printed" << endl;
	// Clean up.
	plot_file.close();
	mat_delete(phi);
	mat_delete(S);
	#endif /* PRACTICE_SESSION */	
	
	
	
	// Practice session 3
	#if PRACTICE_SESSION == 3
	cout << "\nPractice session #3" << endl;
	// MC to be completed.
	#endif /* PRACTICE_SESSION */	
	return 0;
}
