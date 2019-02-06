#include <math.h>
#include <stdio.h>
#define JMAX 20

double qtrap(double (*func)(double), double a, double b, double EPS)
{
	double trapzd(double (*func)(double), double a, double b, int n);
	/*void nrerror(char error_text[]);*/
	int j;
	double s,olds,t1,t2;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=trapzd(func,a,b,j);
		if (fabs(s-olds) < EPS*fabs(olds)) return s;
		olds=s;
	}
	return s;
}
#undef JMAX
