#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../mclib/mclib.h"

#define PRACTICE_SESSION 3

#if PRACTICE_SESSION == 1

#define NX 32
#define NY NX
#define IBEG 1
#define IEND (NX - 2)
#define JBEG 1
#define JEND (NY - 2)

static double const global_S = 0.0;

double solution(double x, double y) {
	return exp(-M_PI * x) * sin(-M_PI * y) + global_S / 4.0 * (x*x + y*y);
}

// Only Dirichlet boundary conditions can be assigned only the first time, Neumann boundary conditions must be applied at each iteration.
void set_boundary_conditions(double ** phi, double x[], double y[]) {
	for (int i = 0; i < NX; i++) {
		phi[i][0] = solution(x[i], y[0]);
		phi[i][NY-1] = solution(x[i], y[NY-1]);
	}
	for (int j = 1; j < NY - 1; j++) {
		phi[0][j] = solution(x[0], y[j]);
		phi[NX-1][j] = solution(x[NX-1], y[j]);
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

#endif /* PRACTICE_SESSION == 1 */
#if PRACTICE_SESSION == 2

#define NX 128
#define NY NX

static double const global_rho_0 = 1.0;
static double const global_a = 0.1;

double rho(double x, double y) {
	double r, _return;
	r = sqrt(x*x + y*y);
	if (r <= global_a) {
		_return = global_rho_0;
	} else {
		_return = 0.0;
	}
	return _return;
}

double solution(double x, double y) {
	double r, _return;
	r = sqrt(x*x + y*y);
	if (r <= global_a) {
		_return = - global_rho_0 * r*r / 4.0;
	} else {
		_return = - global_rho_0 * global_a*global_a / 2.0 * (log(r / global_a) + 0.5);
	}
	return _return;
}

void set_boundary_conditions(double ** phi, double x[], double y[]) {
	for (int i = 0; i < NX; i++) {
		phi[i][0] = solution(x[i], y[0]);
		phi[i][NY-1] = solution(x[i], y[NY-1]);
		// If we know something more on the solution, like its values at the source, we can reduce the number of iterations.
		/*for (int j = 1; j <= NY - 2; j++) {
			// In this case is more efficient to put all the conditions of the y axis together.
			if (i == 0 || i == NX - 1) {
				phi[0][j] = solution(x[0], y[j]);
				phi[NX-1][j] = solution(x[NX-1], y[j]);
			} else if (sqrt(x[i]*x[i] + y[j]*y[j]) <= global_a) {
				phi[i][j] = solution(x[i], y[j]);
			} 
		}*/
	}
	// In case we know the conditions only at the boundaries of the grid, we do not need a nested loop.
	for (int j = 1; j <= NY - 2; j++) {
		phi[0][j] = solution(x[0], y[j]);
		phi[NX-1][j] = solution(x[NX-1], y[j]);
	}
}

double residual(double ** phi, double ** S, double dx, double dy) {
	double epsilon, delta_x, delta_y;
	epsilon = 0.0;
	for (int i = 1; i <= NX - 2; i++) {
		for (int j = 1; j <= NY - 2; j++) {
			delta_x = phi[i+1][j] - 2.0 * phi[i][j] + phi[i-1][j];
			delta_y = phi[i][j+1] - 2.0 * phi[i][j] + phi[i][j-1];
			epsilon += fabs(delta_x + delta_y - dx*dy * S[i][j]);
		}
	}
	return epsilon;
}

#endif /* PRACTICE_SESSION == 2 */
#if PRACTICE_SESSION == 3

#define NX 129
#define NY 65

void set_boundary_conditions(double ** phi, double x[], double y[]) {
	for (int i = 0; i < NX; i++) {
		// Strictly this condition is not necessary because conditions apply to the entire domain of x anyway.
		if (x[i] >= 0 && x[i] <= 2.0) {
			phi[i][0] = 0.0;
			phi[i][NY-1] = 2.0 - x[i];
		}
	}
}

double residual(double ** phi, double ** phi_prev, double dx, double dy) {
	double epsilon = 0.0;
	for (int i = 1; i <= NX - 2; i++) {
		for (int j = 1; j <= NY - 2; j++) {
			epsilon += fabs(phi[i][j] - phi_prev[i][j]) * dx * dy;
		}
	}
	return epsilon;
}

#endif /* PRACTICE_SESSION == 3 */

int main() {
	// Setup.
	using namespace std;
	cout << fixed << setprecision(N_PRECISION);
	double const tolerance = 1e-7;
	int c;
	double x_min, x_max, y_min, y_max, dx, dy, epsilon, omega;
	double ** phi, ** S, x[NX], y[NY];
	ofstream plot_file;
	
	#if PRACTICE_SESSION == 1
	
	cout << "Practice session #1" << endl;
	cout << "S = " << global_S << endl;
	x_min = 0.0;
	x_max = 1.0;
	y_min = 0.0;
	y_max = 1.0;
	dx = (x_max - x_min) / (NX - 1);
	dy = (y_max - y_min) / (NY - 1);
	for (int i = 0; i < NX; i++) {
		x[i] = x_min + i * dx;
	}
	for (int j = 0; j < NY; j++) {
		y[j] = y_min + j * dy;
	}
	
	// Solve with Jacobi method.
	cout << "Jacobi";
	phi = mat_new(NX, NY);
	S = mat_new(NX, NY);
	// Assign internal points and assign boundary conditions.
	mat_zero(phi, NX, NY);
	set_boundary_conditions(phi, x, y);
	mat_constant(S, global_S, NX, NY);
	// Find solution.
	c = 0;
	epsilon = residual(phi, S, dx, dx);
	while (epsilon > tolerance) {
		++c;
		jacobi(phi, S, dx, NX, NY);
		epsilon = residual(phi, S, dx, dx);
	}
	--c;
	cout << ", solved in " << c << " iterations";
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
	cout << ", file printed" << endl;
	// Clean up.
	plot_file.close();
	mat_delete(phi);
	mat_delete(S);
	
	// Solve with Gauss-Seidel method.
	cout << "Gauss-Seidel";
	phi = mat_new(NX, NY);
	S = mat_new(NX, NY);
	// Assign internal points and assign boundary conditions.
	mat_zero(phi, NX, NY);
	set_boundary_conditions(phi, x, y);
	mat_constant(S, global_S, NX, NY);
	// Find solution.
	c = 0;
	epsilon = residual(phi, S, dx, dx);
	while (epsilon > tolerance) {
		++c;
		gauss_seidel(phi, S, dx, NX, NY);
		epsilon = residual(phi, S, dx, dx);
	}
	--c;
	cout << ", solved in " << c << " iterations";
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
	cout << ", file printed" << endl;
	// Clean up.
	plot_file.close();
	mat_delete(phi);
	mat_delete(S);
	
	// Solve with successive over relaxation method.
	cout << "SOR";
	omega = 2.0 / (1.0 + M_PI / NX);
	phi = mat_new(NX, NY);
	S = mat_new(NX, NY);
	// Assign internal points and assign boundary conditions.
	mat_zero(phi, NX, NY);
	set_boundary_conditions(phi, x, y);
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
	cout << ", solved in " << c << " iterations";
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
	cout << ", file printed" << endl;
	// Clean up.
	plot_file.close();
	mat_delete(phi);
	mat_delete(S);
	
	#endif /* PRACTICE_SESSION == 1 */
	#if PRACTICE_SESSION == 2

	cout << "Practice session #2" << endl;
	// Grid preparation.
	x_min = -1.0;
	x_max = 1.0;
	y_min = -1.0;
	y_max = 1.0;
	dx = (x_max - x_min) / (NX - 1);
	dy = (y_max - y_min) / (NY - 1);
	for (int i = 0; i < NX; i++) {
		x[i] = x_min + i * dx;
	}
	for (int j = 0; j < NY; j++) {
		y[j] = y_min + j * dy;
	}

	// Solve with Jacobi method.
	cout << "Jacobi";
	phi = mat_new(NX, NY);
	S = mat_new(NX, NY);
	// Assign initial values and boundary conditions.
	mat_zero(phi, NX, NY);
	set_boundary_conditions(phi, x, y);
	for (int i = 0; i < NX; i++) {
		for (int j = 0; j < NY; j++) {
			S[i][j] = -rho(x[i], y[j]);
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
	cout << ", solved in " << c << " iterations";
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
	cout << ", file printed" << endl;
	// Clean up.
	plot_file.close();
	mat_delete(phi);
	mat_delete(S);

	// Solve with Gauss-Seidel method.
	cout << "Gauss-Seidel";
	phi = mat_new(NX, NY);
	S = mat_new(NX, NY);
	// Assign initial values and boundary conditions.
	mat_zero(phi, NX, NY);
	set_boundary_conditions(phi, x, y);
	for (int i = 0; i < NX; i++) {
		for (int j = 0; j < NY; j++) {
			S[i][j] = -rho(x[i], y[j]);
		}
	}
	// Find solution.
	c = 0;
	epsilon = residual(phi, S, dx, dx);
	while (epsilon > tolerance) {
		++c;
		gauss_seidel(phi, S, dx, NX, NY);
		epsilon = residual(phi, S, dx, dx);
	}
	--c;
	cout << ", solved in " << c << " iterations";
	// Print results on file.
	plot_file.open("elliptic_2_gauss_seidel.dat");
	plot_file << fixed << setprecision(N_PRECISION);
	plot_file << "x y phi(x,y)" << endl;
	for (int j = 0; j < NY; j++) {
		for (int i = 0; i < NX; i++) {
			plot_file << x[i] << ' ' << y[j] << ' ' << phi[i][j] << endl;
		}
		plot_file << endl;
	}
	cout << ", file printed" << endl;
	// Clean up.
	plot_file.close();
	mat_delete(phi);
	mat_delete(S);

	// Solve with successive over relaxation method.
	cout << "SOR";
	omega = 2.0 / (1.0 + M_PI / NX);
	phi = mat_new(NX, NY);
	S = mat_new(NX, NY);
	// Assign internal points and assign boundary conditions.
	mat_zero(phi, NX, NY);
	set_boundary_conditions(phi, x, y);
	for (int i = 0; i < NX; i++) {
		for (int j = 0; j < NY; j++) {
			S[i][j] = -rho(x[i], y[j]);
		}
	}
	// Find solution.
	c = 0;
	epsilon = residual(phi, S, dx, dx);
	while (epsilon > tolerance) {
		++c;
		successive_over_relaxation(phi, S, dx, omega, NX, NY);
		epsilon = residual(phi, S, dx, dx);
	}
	--c;
	cout << ", solved in " << c << " iterations";
	// Print results on file.
	plot_file.open("elliptic_2_SOR.dat");
	plot_file << fixed << setprecision(N_PRECISION);
	plot_file << "x y phi(x,y)" << endl;
	for (int j = 0; j < NY; j++) {
		for (int i = 0; i < NX; i++) {
			plot_file << x[i] << ' ' << y[j] << ' ' << phi[i][j] << endl;
		}
		plot_file << endl;
	}
	cout << ", file printed" << endl;
	// Clean up.
	plot_file.close();
	mat_delete(phi);
	mat_delete(S);

	#endif /* PRACTICE_SESSION == 2 */
	#if PRACTICE_SESSION == 3
	
	cout << "Practice session #3" << endl;
	// Grid preparation.
	x_min = 0.0;
	x_max = 2.0;
	y_min = 0.0;
	y_max = 1.0;
	dx = (x_max - x_min) / (NX - 1);
	dy = (y_max - y_min) / (NY - 1);
	for (int i = 0; i < NX; i++) {
		x[i] = x_min + i * dx;
	}
	for (int j = 0; j < NY; j++) {
		y[j] = y_min + j * dy;
	}
	
	// Solve with successive over relaxation method.
	cout << "SOR";
	double sigma, phi_prev;
	omega = 2.0 / (1.0 + M_PI / NX);
	phi = mat_new(NX, NY);
	S = mat_new(NX, NY);
	// Assign initial values and boundary conditions at the grid boundaries.
	mat_zero(phi, NX, NY);
	set_boundary_conditions(phi, x, y);
	mat_zero(S, NX, NY);
	// Find solution.
	c = 0;
	do {
		++c;
		epsilon = 0.0;
		// Successive over relaxation internal loop is modified to implement Von Neumann boundary conditions and the new formula for the residual.
		for (int i = 1; i <= NX - 2; i++) {
			for (int j = 1; j <= NY - 2; j++) {
				phi_prev = phi[i][j];
				// Notice that these formulas are valid for dx = dy.
				if (i == 1) {
					sigma = 0.0;
					phi[i][j] = ( (1.0 - omega) * phi[i][j] + omega / 4.0 * (-dx * sigma + phi[i+1][j] + phi[i][j-1] + phi[i][j+1] - dx*dx * S[i][j]) ) / (1.0 - omega / 4.0);
				} else if (i == NX - 2) {
					sigma = 3.0;
					phi[i][j] = ( (1.0 - omega) * phi[i][j] + omega / 4.0 * (phi[i-1][j] + dx * sigma + phi[i][j-1] + phi[i][j+1] - dx*dx * S[i][j]) ) / (1.0 - omega / 4.0);
				} else {
					phi[i][j] = (1.0 - omega) * phi[i][j] + omega / 4.0 * (phi[i-1][j] + phi[i+1][j] + phi[i][j-1] + phi[i][j+1] - dx*dx * S[i][j]);
				}
				// Evaluate residual.
				epsilon += fabs(phi[i][j] - phi_prev) * dx * dy;
			}
		}
	} while (epsilon > tolerance);
	--c;
	cout << ", solved in " << c << " iterations";
	// Print results on file.
	plot_file.open("elliptic_3_SOR.dat");
	plot_file << fixed << setprecision(N_PRECISION);
	plot_file << "x y phi(x,y)" << endl;
	for (int j = 0; j < NY; j++) {
		for (int i = 0; i < NX; i++) {
			plot_file << x[i] << ' ' << y[j] << ' ' << phi[i][j] << endl;
		}
		plot_file << endl;
	}
	cout << ", file printed" << endl;
	// Clean up.
	plot_file.close();
	mat_delete(phi);
	mat_delete(S);
	
	#endif /* PRACTICE_SESSION == 3 */

	return 0;
}
