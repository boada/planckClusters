#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "header.h"


void output_halo_mass_function()
{
  int k,nr=100;
  double x,dlogm,m,mvir,cdelta,mmin=1e8,mmax=1e16,delta_vir;
  FILE *fp;
  char aa[100];

  fprintf(stderr,"\n\nCALCULATING HALO MASS FUNCTION.\n");
  fprintf(stderr,    "-------------------------------\n\n");

  sprintf(aa,"%s.dndM",Task.root_filename);
  fp = fopen(aa,"w");

  dlogm = (log(mmax) - log(mmin))/(nr-1);

  for(k=0;k<nr;++k)
    {
      m = exp(k*dlogm)*mmin;
      x = halo_mass_function(m);
      fprintf(fp,"%e %e\n",m,x);
    }
  fclose(fp);
}
