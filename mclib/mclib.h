#include <cmath>
#include <iostream>



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

