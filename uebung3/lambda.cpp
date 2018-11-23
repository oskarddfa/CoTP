#define _USE_MATH_DEFINES

#include <math.h>

double lambda(double a, double A[], double b, double B[])
	{
	return pow((M_PI/(a+b)),1.5)*exp(-a*b/(a+b)*((A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])+(A[2]-B[2])*(A[2]-B[2]));
	}
