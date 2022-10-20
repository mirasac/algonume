#include <cmath>
#include <iostream>

#define DEBUG



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
void gaussianpoints(int Ng, double * t, double * w);

/*
Check if intermediate value theorem is applicable for function f in the interval [a, b], in other words check if the function, if continuous, may have a zero in the interval.

Return
------
Value 1 if the theorem is valid for the given parameters, 0 otherwise.
*/
char check_intermediate_value(double (*f)(double x), double a, double b);



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
