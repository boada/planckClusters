#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "header.h"

/* This calculates the linear growthfactor at redshift z.
 * Currently, the default is a flat LCDM universe, but can easily
 * be changed.
 */
double x03_g1,
  xi3_g1;
double func_D0(double);
double func_Di(double);

double growthfactor(double z)
{
  int i;
  double zp1,x03,xi3,bi,b0,yy,sqx,sqx1,dbdx,x0,xi,lambda_i,htemp,omega_L,omega_i,hubble_i,
    astart,fac,redshift;

  redshift=z;
  astart=1.0/(z+1);
  zp1=redshift+1;
  omega_L=1-OMEGA_M;

  hubble_i = sqrt(OMEGA_M/(astart*astart*astart) + 
		  (1.0-OMEGA_M-omega_L)/(astart*astart) +omega_L);

  if(omega_L>0)
    {
      lambda_i=omega_L/(omega_L+(1.0-OMEGA_M-omega_L)*zp1*zp1+OMEGA_M*pow(zp1,3.0));
      omega_i=OMEGA_M*pow(zp1,3.0)*lambda_i/omega_L;
    }
  else
    {
      lambda_i=0;
      omega_i=OMEGA_M*zp1/(1.0+OMEGA_M*redshift);
    }

  fac=astart;

  if((OMEGA_M < 0.99) && (omega_L > 0.001))
    {
      x03_g1=x03=1.0/OMEGA_M-1.0;
      xi3_g1=xi3=1.0/omega_i-1.0;
      b0 = qromo(func_D0,0.0,1.0,midpnt);
      bi = qromo(func_Di,0.0,1.0,midpnt);
      b0=b0*sqrt(x03+1.0);
      bi=bi*sqrt(xi3+1.0)*astart;
      fac=bi/b0;
    }


  if((OMEGA_M < 0.99) && (omega_L == 0))
    {
      x0=1.0/OMEGA_M-1.0;
      xi=x0*astart;
      b0 = 1.+3./x0+3.*sqrt(1+x0)*log(sqrt(1.+x0)-sqrt(x0))/pow(x0,1.5);
      bi = 1.+3./xi+3.*sqrt(1+xi)*log(sqrt(1.+xi)-sqrt(xi))/pow(xi,1.5);
      fac = bi/b0;
    }
  return(fac);
}

double func_Di(double y)
{
  return(pow(1.0+xi3_g1*pow(y,1.2),-1.5));
}

double func_D0(double y)
{
  return(pow(1.0+x03_g1*pow(y,1.2),-1.5));
}
