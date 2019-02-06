#include <math.h>
#include <stdio.h>
#include "header.h"
#define EPS 1.0e-6
#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5

double qromo(double (*func)(double), double a, double b,
	double (*choose)(double(*)(double), double, double, int))
{
	void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
	void nrerror(char error_text[]);
	int j;
	double ss,dss,h[JMAXP+1],s[JMAXP+1];

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=h[j]/9.0;
	}
	printf("ERROR: Too many steps in routing qromo\n");
	ERROR_FLAG=1;
	return ss;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K
