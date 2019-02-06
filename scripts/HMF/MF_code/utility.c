#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "header.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

#define NR_END 1
#define FREE_ARG char*

/* returns the number of cpu-ticks in seconds that
 * have elapsed. (or the wall-clock time)
 */
double second(void)
{
  return ((double)((unsigned int)clock()))/CLOCKS_PER_SEC;

  /* note: on AIX and presumably many other 32bit systems, 
   * clock() has only a resolution of 10ms=0.01sec 
   */ 
}


/* returns the time difference between two measurements 
 * obtained with second(). The routine takes care of the 
 * possible overflow of the tick counter on 32bit systems.
 */
double timediff(double t0,double t1)
{
  double dt;
  
  dt=t1-t0;

  if(dt<0)  /* overflow has occured */
    {
      dt=t1 + pow(2,32)/CLOCKS_PER_SEC - t0;
    }

  return dt;
}


void endrun(char *instring)
{
  fprintf(stderr,"endrun> %s\n",instring);
  fflush(stderr);
  /*
  fprintf(stdout,"endrun> %s\n",instring);
  fflush(stdout);
  */

#ifdef PARALLEL
  MPI_Abort(MPI_COMM_WORLD, 0);
#endif
  exit(0);


}

/* This takes a file and reads the number of lines in it,
 * rewinds the file and returns the lines.
 */
int filesize(FILE *fp)
{
  int i=-1;
  char a[1000];

  while(!feof(fp))
    {
      i++;
      fgets(a,1000,fp);
    }
  rewind(fp);
  return(i);
}

/* This opens a file and has an error trap
 * if the file does not exist.
 */
FILE *openfile(char *ff)
{
  FILE *fp;
  if(!(fp=fopen(ff,"r")))
    {
      fprintf(stderr,"ERROR opening [%s]\n",ff);
      exit(0);
    }
  return(fp);
}


/* This is for allocating a 3-dimensional array of doubles, adapted from
 * the numerical recipes f3tensor routine.
 */
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

/* Free the memory allocated by d3tensor (see above)
 */
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}
