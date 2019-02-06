#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
	void nrerror(char error_text[]);
	int klo,khi,k,i;
	double h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if(h == 0.0) {
	  printf("BADXA %d %f\n",n,x);
	  fflush(stdout);
	}
	if (h == 0.0) nrerror("Bad xa input to routine splint");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
	if(isnan(*y))
	  {
	    /*
	    fprintf(stdout,"NAN %e %e %d %e %e %e %e %e %e\n",a,b,n,x,ya[klo],ya[khi],
		    y2a[khi],y2a[klo],h);
	    for(i=1;i<=n;++i)
	      printf("NAN %d %e %e %e %e\n",i,xa[i],ya[i],y2a[i],x);
	    */
	    *y=0;
	  }
}
