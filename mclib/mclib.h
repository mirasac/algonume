#ifndef __MCLIB_H
#define __MCLIB_H

#include <cmath>
#include <iostream>

#define N_PRECISION 12
#define FLG_DEBUG 0



////////// Utilities functions //////////



/*
Set values a and b in ascending order.
*/
void orderinterval(double * a, double * b);

/*
Get the right gaussian points and weigths given their number. Gaussian points x are given in ascending order and the respective weights are ordered accordingly.

Parameters
----------
Ng : int
	Number of gaussian points.
x : double *
	Array of length Ng filled with the gaussian points.
w : double *
 	Array of length Ng filled with weights associated to each gaussian point.
*/
void gaussianpoints(int Ng, double t[], double w[]);

/*
Check if intermediate value theorem is applicable for function f in the interval [a, b], in other words check if the function, if continuous, may have a zero in the interval.

Return
------
Value 1 if the theorem is valid for the given parameters, 0 otherwise.
*/
char check_intermediate_value(double (*f)(double x), double a, double b);

/*
Evaluate the polynomial of order n defined by coefficients c in point x.

Parameters
----------
n : int
	Order of the polynomial.
c : double[]
	Array of coefficients, coefficient in position i is multiplied to x^i.
x : double
	Point where the polynomial is evaluated.
*/
double polynomial(int n, double c[], double x);

/*
Set array c1 as the array of coefficients of the derivative of the polynomial of coefficients c.
*/
void polynomial_derivative(int n, double c[], double c1[]);

/*
Evaluate Legendre polynomial of order n in point x.
*/
double polynomial_legendre(int n, double x);

// MC does not work because I do not know how to use correctly the preprocessor, elaborate.
/*
#define POW2(x) ( (x)*(x) )

#define POW3(x) ( (x)*(x)*(x) )

#define POW4(x) ( (x)*(x)*(x)*(x) )

#define POW5(x) ( (x)*(x)*(x)*(x)*(x) )

#define POW6(x) ( (x)*(x)*(x)*(x)*(x)*(x) )

#define POW7(x) ( (x)*(x)*(x)*(x)*(x)*(x)*(x) )

#define POW8(x) ( (x)*(x)*(x)*(x)*(x)*(x)*(x)*(x) )

#define POW9(x) ( (x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x) )

#define POW(x, n) POW n(x)
*/


////////// Quadrature functions //////////



/*
Evaluate the integral of function f in the real interval [a, b] divided in N sub intervals using the rectangular quadrature.
*/
double rectangularquad(double (*f)(double x), double a, double b, int N);

/*
Evaluate the integral of function f in the real interval [a, b] divided in N sub intervals using the midpoint quadrature.
*/
double midpointquad(double (*f)(double x), double a, double b, int N);

/*
Evaluate the integral of function f in the real interval [a, b] divided in N sub intervals using the trapezioidal quadrature.
*/
double trapezioidalquad(double (*f)(double x), double a, double b, int N);

/*
Evaluate the integral of function f in the real interval [a, b] divided in N sub intervals using the Simpson rule.
*/
double simpsonquad(double (*f)(double x), double a, double b, int N);

/*
Evaluate the integral of function f in the real interval [a, b] divided in N sub intervals using Gauss-Legendre quadrature.
*/
double gaussquad(double (*f)(double x), double a, double b, int N, int Ng);

/*
Evaluate the integral of function f in the real interval [a, b] x [a, b] each divided in N sub intervals using Gauss-Legendre quadrature.
*/
double multiquad(double (*f)(double x, double y), double a, double b, int N, int Ng);



////////// Root finders //////////



/*
Find a zero of function f in the interval [a, b] within the specified tollerance using the bisection method.

Return
------
The first point x where the function assumes a value lower than the specified tollerance.
*/
double bisection(double (*f)(double x), double a, double b, double tollerance);

/*
Find a zero of function f in the interval [a, b] within the specified tollerance using the false position method.

Return
------
The first point x where the function assumes a value lower than the specified tollerance.
*/
double falseposition(double (*f)(double x), double a, double b, double tollerance);

/*
Find a zero of function f in the interval [a, b] within the specified tollerance using the secant method.

Return
------
The first point x where the function assumes a value lower than the specified tollerance.
*/
double secant(double (*f)(double x), double a, double b, double tollerance);

/*
Find a zero of function f in the interval [a, b] within the specified tollerance using the Newton-Raphson method.

Return
------
The first point x where the function assumes a value lower than the specified tollerance.
*/
double newtonraphson(double (*f)(double x), double (*f1)(double x), double a, double b, double tollerance);

/*
Variant of function `newtonraphson` to be used when function f is a polynomial.
*/
double newtonraphson_poly(double (*p)(int n, double c[], double x), int n, double c[], double a, double b, double tollerance);



////////// Bracketing //////////



/*
Execute the bracketing of interval [a, b] by dividing it in N subintervals and checking for each subinterval if function f changes sign. Boundaries of intervals which satisfy this condition are saved in arrays x_L (lower boundaries) and x_R (upper boundaries) in ascending order.

Return
------
Number of intervals which may have at least a zero of the function.
*/
int bracketing(double (*f)(double x), double a, double b, int N, double x_L[], double x_R[]);



////////// Differentation //////////



/*
Evaluate forward difference of f in point x with grid step h.
*/
double forward_difference(double (*f)(double), double h, double x);

/*
Evaluate backward difference of f in point x with grid step h.
*/
double backward_difference(double (*f)(double), double h, double x);

/*
Evaluate central difference of f in point x with grid step h.
*/
double central_difference(double (*f)(double), double h, double x);



////////// ODE integration //////////



/*
Perform a single step of the Euler method for ODEs resolution. First order ODEs are accepted, but an ODE of any order can be treated if it is rewrited as a system of first order ODEs.

Parameters
----------
t : double
	Time value considered during the step execution.
dt : double
	Time step.
Y : double
	Array of numerical solutions at time t + dt.
rhs : void
	Function representing the RHS terms of the equations. It operates simultaneously on all equations, using t as time, Y_0 as array of initial values and R as array of solutions of the operation.
n_eq : int
	Number of first order ODEs treated.
*/
void eulerstep(double t, double dt, double Y[], void (*rhs)(double t, double Y_0[], double R[]), int n_eq);

#endif /* __MCLIB_H */
