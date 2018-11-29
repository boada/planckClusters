#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "header.h"

/* This function is taken from Volker Springel's GADGET code and
 * modified for the parameters used for the HOD.x code.
 */

/*
 *  This function parses the parameterfile in a simple way.
 *  Each paramater is defined by a keyword (`tag'), and can be
 *  either of type douple, int, or character string.
 *  The routine makes sure that each parameter appears 
 *  exactly once in the parameterfile.
 */
void read_parameter_file(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define CHAR 4
#define LONG 4
#define MAXTAGS 300

  FILE *fd,*fdout;

  char buf[200],buf1[200],buf2[200],buf3[200],tempchar;
  int  i,j,nt,ii,nn,ctemp;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][200];
  int  errorFlag=0;
  int IDUM_MCMC_TEMP=-555;

  nt=0;
  
  strcpy(tag[nt],"REDSHIFT");
  addr[nt]=&REDSHIFT;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"GAMMA");
  addr[nt]=&GAMMA;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"OMEGA_M");
  addr[nt]=&OMEGA_M;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"OMEGA_B");
  addr[nt]=&OMEGA_B;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"SIGMA_8");
  addr[nt]=&SIGMA_8;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"RHO_CRIT");
  addr[nt]=&RHO_CRIT;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"HUBBLE");
  addr[nt]=&HUBBLE;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"ITRANS");
  addr[nt]=&ITRANS;
  id[nt++]=INT;
    
  strcpy(tag[nt],"SPECTRAL_INDX");
  addr[nt]=&SPECTRAL_INDX;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"DELTA_CRIT");
  addr[nt]=&DELTA_CRIT;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"DELTA_HALO");
  addr[nt]=&DELTA_HALO;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"TF_file"); 
  addr[nt]=Files.TF_file;
  id[nt++]=STRING;

 strcpy(tag[nt],"root_filename");
  addr[nt]=&Task.root_filename;
  id[nt++]=STRING;
   
  if((fd=fopen(fname,"r")))
    {
      nn=filesize(fd);
      sprintf(buf,"%s","hod-usedvalues");
      if(!(fdout=fopen(buf,"w")))
	{
	  fprintf(stdout,"error opening file '%s' \n",buf);
	  errorFlag=1; 
	}
      else
	{
	  /*while(!feof(fd))*/
	  for(ii=1;ii<=nn;++ii)
	    {
	      fgets(buf,200,fd);
	      if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
		continue;
	      
	      if(buf1[0]=='%')
		continue;
	      
	      for(i=0,j=-1;i<nt;i++)
		if(strcmp(buf1,tag[i])==0)
		  {
		    j=i;
		    tag[i][0]=0;
		    break;
		  }
	      
	      if(j>=0)
		{
		  switch(id[j])
		    {
		    case DOUBLE:
		      *((double*)addr[j])=atof(buf2); 
		      fprintf(fdout,"%-35s%g\n",buf1,*((double*)addr[j]));
		      break;
		    case STRING:
		      strcpy(addr[j],buf2);
		      fprintf(fdout,"%-35s%s\n",buf1,buf2);
		      break;
		    case INT:
		      *((int*)addr[j])=atoi(buf2);
		      fprintf(fdout,"%-35s%d\n",buf1,*((int*)addr[j]));
		      break;
		    case CHAR:
		      *((char*)addr[j])=buf2[0];
		      fprintf(fdout,"%-35s%c\n",buf1,*((int*)addr[j]));
		      break;
		    }
		}
	      else
		{
		  fprintf(stderr,"Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			  fname,buf1);
		  errorFlag=1;
		}
	    }
	}
      fclose(fd);
      fclose(fdout);

    }
  else
    {
      fprintf(stderr,"Parameter file %s not found.\n", fname);
      exit(1);
    }



  for(i=0;i<nt;i++)
    {      
      if(!strcmp(tag[i],"GAMMA"))
	{
	  if(ITRANS==4) {
	    fprintf(stderr,"No value of GAMMA given for ITRANS=4\n");
	    errorFlag=1;
	  }
	  continue;
	}
      if(!strcmp(tag[i],"TF_file"))
	{
	  if(ITRANS==11) {
	    sprintf(Files.TF_file,"CMBFAST_trans.dat");
	    fprintf(stderr,"No transfer function file, using [%s]\n",Files.TF_file);
	  }
	  continue;
	}

      if(*tag[i])
	{
	  fprintf(stderr,"Error. I miss a value for tag '%s' in parameter file '%s'.\n",
		  tag[i],fname);
	  errorFlag=1;
	}
    }

  if(errorFlag)
    endrun("error in input_params ");

  /* If REDSHIFT>0, calculate the amplitude of fluctuations at the desired z
   */
  if(REDSHIFT>0)
    {
      SIGMA_8 = SIGMA_8*growthfactor(REDSHIFT);
      fprintf(stderr,"SIGMA_8 @ z=%f is %f\n",REDSHIFT,SIGMA_8);
    }

  /* Other initialization stuff.
   */
  MSTAR=mstar();

#undef DOUBLE 
#undef STRING 
#undef INT 
#undef MAXTAGS
}
 

