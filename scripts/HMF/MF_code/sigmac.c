#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

/*
 *    sigmac calls qromo to evaluate the integral

         int_0^infty P_s(k) 4 pi k^2 dk / (2pi)^3

      where P_s(k) is the power spectrum smoothed through a gaussian
      window of radius rgaus and a top hat window of radius rth.
      The power spectrum is a power law of index xindx, multiplied by
      the square of the CDM transfer function if itrans = 1.

      If the smoothing radius [rad] is positive, a top-hat window
      function is used. 
      If [rad<0], a gaussian window function is used.

*/


double rad1;
double func_sigmac(double xk);

double sigmac_interp(double m)
{
  static int flag=0,prev_cosmo=0;
  static double *x,*y,*y2, pnorm;
  int i,n=100;
  double dlogm,max=5.e16,min=1.0e9,a,b,m1,m2,dm1,xi,power,rm,sig,b1,b2,mass;

  if(!flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      if(!ThisTask && OUTPUT)
	fprintf(stdout,"RESET: resetting bias for %f %f\n",OMEGA_M,SIGMA_8);
      if(!flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      flag=1;
      dlogm = (log(max) - log(min))/(n-1);
      pnorm=SIGMA_8/sigmac(8.0);
      for(i=1;i<=n;++i)
	{
	  m = exp((i-1)*dlogm)*min;
	  rm=pow(3.0*m/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
	  sig=log(pnorm*sigmac(rm));
	  x[i] = log(m);
	  y[i] = sig;
	}
      spline(x,y,n,2.0E+30,2.0E+30,y2);
      prev_cosmo=RESET_COSMOLOGY;

    }
  m=log(m);
  splint(x,y,y2,n,m,&a);
  return exp(a);
}


double sigmac2(double rad)
{
  int i;
  double xk1,s1=1,s2=0,sigma0;

  rad1=rad;
  xk1 = 1./(2.*PI*mabs(rad)); 
  s1=qromo(func_sigmac,0.0,xk1,midpnt);
  for(s2=0,i=1;i<=10;++i)
    {
      //s2 += qromo(func_sigmac,xk1,16*xk1,midpnt);
      s2 += qtrap(func_sigmac,xk1,16*xk1,1.0E-5);
      xk1 = 16*xk1;
      printf("%d %e\n",i,s2);
    }
  printf("boo %e %e\n",s1,s2);
  sigma0=sqrt((s1+s2)*(4*PI)/pow(2*PI,3.0));
  return sigma0; 
  
}
double sigmac(double rad)
{
  double xk1,s1=1,s2=0,sigma0;

  rad1=rad;
  xk1 = 1./(2.*PI*mabs(rad)); 
  s1=qromo(func_sigmac,0.0,xk1,midpnt);
  s2=qromo(func_sigmac,xk1,1.e20,midinf);
  sigma0=sqrt((s1+s2)*(4*PI)/pow(2*PI,3.0));
  return sigma0; 
  
}

double func_sigmac(double xk)
{
  double xkr,w,psp;

  if(rad1>0)
    {
      xkr = xk*rad1;
      w = 3*(sin(xkr)-xkr*cos(xkr))/(xkr*xkr*xkr);
      w = w*w;
    }
  else
    {
      xkr = -rad1*xk;
      w = exp(-xkr*xkr);
    }

  if(ITRANS>0)
    psp=pow(xk,SPECTRAL_INDX)*pow(transfnc(xk),2.0);
  else
    psp=pow(xk,SPECTRAL_INDX);
  psp=psp*w*xk*xk;
  return(psp);
}

