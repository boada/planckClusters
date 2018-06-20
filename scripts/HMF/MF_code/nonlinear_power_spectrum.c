#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

/* This is my implementation of the Halo Model fitting function
 * described in Appendix C of Smith, Peacock, et al. 2003 MNRAS 341, 1311
 *
 * Comparison to Smith et al. program halofit.f output is successful.
 */

/* This returns the linear power 
 *   Delta = 4pi*k^3 P(k)/(2pi)^3
 */
double linear_power_spectrum(double xk)
{
  static double *kk,*pknl,*y2,pnorm=-1,ahi,bhi;
  static int flag=1,nk=1000,prev_cosmology=0;
  double a,psp,x1[4],y1[4];
  int i;

  if(pnorm<0 || prev_cosmology!=RESET_COSMOLOGY)
    {
      pnorm=SIGMA_8/sigmac(8.0);  
      pnorm*=pnorm;
      prev_cosmology=RESET_COSMOLOGY;
    }
  if(ITRANS>0)
    psp=pow(xk,SPECTRAL_INDX)*pow(transfnc(xk),2.);
  else
    psp=pow(xk,SPECTRAL_INDX);
  //  printf("BOO %e %e\n",xk,psp*pnorm);
  psp=psp*pnorm*xk*xk*xk/(2*PI*PI);
  return(psp);
}

