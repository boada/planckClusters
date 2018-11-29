#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "header.h"

/* This routine reads in a transfer function from a file. 
 * - col 1 = k [h/Mpc]
 * - col 2 = T(k)
 * other columns not used.
 *
 * The TF is considered to be un-normalized, and will normalize all values 
 * by the entry in the first line.
 *
 * The TF is stored in arrays an spline interpolation is used at all k values.
 * The interpolation is in log(k)/log(TF) to preserve the power-law dependence
 * of T(k) on k outside the k-range of the file.
 */

double transfunc_file(double xk)
{
  static double *x,*y,*y2, aa, bb;
  static int flag=1,n;
  int i;
  double t,x0;
  FILE *fp;
  char a[1000];
  float x1,x2;

  if(flag)
    {
      flag=0;

      fp=openfile(Files.TF_file);
      n=filesize(fp);

      x=dvector(1,n);
      y=dvector(1,n);
      y2=dvector(1,n);
      for(i=1;i<=n;++i)
	{
	  fscanf(fp,"%f %f",&x1,&x2);
	  x[i]=x1;
	  y[i]=x2;
	  if(i==1){ x0=y[i]; }
	  fgets(a,1000,fp);
	  y[i]/=x0;
	  x[i]=log(x[i]);
	  y[i]=log(y[i]);
	  //printf("BOO %f %f\n",x[i],y[i]);
	}
      fclose(fp);
      spline(x,y,n,1.0E+30,1.0E+30,y2);
      bb = (y[n-5] - y[n])/(x[n-5] - x[n]);
      aa = y[n] - bb*x[n];
    }
  xk=log(xk);
  if(xk>x[n]) 
    return(exp(aa+xk*bb));
  splint(x,y,y2,n,xk,&t);
  return(exp(t));

}


