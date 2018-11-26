#define _USE_MATH_DEFINES

#include <math.h>

double lambda(double a, double b, double Ax, double Bx, double Ay=0.0, double Az=0.0, double By=0.0, double Bz=0.0)
	{
	return pow((M_PI/(a+b)),1.5) * exp( - a*b/(a+b) * ((Ax-Bx)*(Ax-Bx) + (Ay-By)*(Ay-By) + (Az-Bz)*(Az-Bz)));
	}
