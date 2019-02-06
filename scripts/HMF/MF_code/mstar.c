#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "header.h"

/* This function calculates M_star (non-linear mass scale) for the given
 * cosmological paramters.
 */

double sigma_Mmdelta_c(double lnM);
double pnorm1;

double mstar()
{
  double sig,lnMmin,lnMmax,M_star;

  sig=sigmac(8.0); 
  pnorm1 = SIGMA_8/sig;
  
  lnMmin=log(1e3);
  lnMmax=log(1e18);
  M_star=zbrent(sigma_Mmdelta_c,lnMmin,lnMmax,1e-5);
  M_star=exp(M_star);
  if(!ThisTask)
    fprintf(stderr,"M_star = %e h^{-1}M_sol\n",M_star);
  return(M_star); 
}

/*** solve for M_* ***/
double sigma_Mmdelta_c(double lnM)
{ 
  double sig,M,rm;

  M=exp(lnM);
  rm=pow(3.0*M/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  sig=pnorm1*sigmac(rm);

  return sig-DELTA_CRIT;
}
