#include "nrutil.h"
#include "stdlib.h"
#include "stdio.h"

/* Function prototypes--> utility files
 */
double second(void);
double timediff(double t0,double t1);
void endrun(char *instring);
FILE *openfile(char *ff);
int filesize(FILE *fp);
void least_squares(double *x, double *y, int n, double *a, double *b);
void check_for_smoothness(double *x, double *y, int n, double r);
void read_parameter_file(char *fname);
void output_parameter_file(char *fname);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh);

/* Function prototypes--> numerical recipes
 */
double qromo(double (*func)(double), double a, double b,
	     double (*choose)(double(*)(double), double, double, int));
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double midpnt(double (*func)(double), double a, double b, int n);
double midinf(double (*funk)(double), double aa, double bb, int n);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
double zbrent(double (*func)(double), double x1,double x2, double tol);
double qtrap(double (*func)(double), double a, double b, double EPS);
void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret,
	    double (*func)(double []));
void amoeba(double **p, double y[], int ndim, double ftol,
	    double (*funk)(double []), int *nfunk);
void gaussj(double **a, int n, double **b, int m);
float gasdev(long *idum);
void jacobi(double **a, int n, double d[], double **v, int *nrot);
float ran1(long *idum);
double ran2(long *idum);
double trapzd(double (*func)(double), double a, double b, int n);
void sort2(unsigned long n, float arr[], int id[]);


/* Function prototypes--> power spectrum routines.
 */
double transfnc(double xk);
double sigmac(double rad);
double nonlinear_sigmac(double rad);
double transfunc_file(double xk);
double nonlinear_power_spectrum(double xk);
double linear_power_spectrum(double xk);
double mstar(void);
double tf_eisenstein_hu(double k);
double cobenorm(double omega_m);
double cobe_prior(double omega_m);
double sigmac_interp(double m);

/* Function protoypes--> outputting matter stats
 */
void output_matter_power_spectrum(void);
void output_matter_correlation_function(void);
void output_matter_variance(void);
void output_halo_concentrations(void);
void output_halo_mass_function(void);


/* Function prototypes--> halo mass function.
 */
double halo_mass_function(double mass);
double dndM_interp(double m);


/* Function prototypes--> halo concentrations.
 */
double growthfactor(double z);
double halo_concentration(double m);
double cvir_model(double mass);
double halo_mass_conversion(double mvir, double *cvir1, double delta_halo);
double halo_mass_conversion2(double mvir, double cvir1, double delta_vir, double delta_halo);
double halo_c200(double m);


/* Definitions
 */
#define PI       3.14159265358979323846
#define TWOPI    (2*PI)
#define THIRD    (1.0/3.0)
#define ROOT2    1.41421356237309504880
#define RT2PI    2.50662827463100050241
#define LN_2     0.6931471805599452
#define ROOT8    2.82842712475
#define WORKBUF  1000
#define LOGE_10  2.30258509
#define BIG_G    4.301e-9


#define mabs(A)  ((A) < 0.0 ? -(A) : (A))
#define cnint(x) ((x-floor(x)) < 0.5 ? floor(x) : ceil(x))
#define muh(x)   fprintf(stdout,"%d\n",x);fflush(stdout)
#define fmuh(x)  fprintf(stderr,"%e\n",x)
#define square(x) (x*x)

/* Global variables
 */
extern double 
  GAMMA,
  HUBBLE,
  SIGMA_8,
  RHO_CRIT,
  SPECTRAL_INDX,
  OMEGA_M,
  OMEGA_B,
  OMEGA_TEMP,
  DELTA_CRIT,
  MSTAR,
  REDSHIFT,
  DELTA_HALO,
  JENKINS_A,
  JENKINS_B,
  JENKINS_C;

extern int RESET_COSMOLOGY,
  ITRANS,
  OUTPUT,
  ERROR_FLAG;

extern int ThisTask;

/* Various input files and flags on whether or not to use them.
 */
extern struct file_parameters {
  char HaloFile[1000];
  char HaloDensityFile[1000];
  char TF_file[100];
  char TwoHaloFile[100];
  int  i_TwoHalo;
  char MassFuncFile[100];
  int i_MassFunc;
  char PDFTable[100];
  int i_PDFTable;
  char PSPFile[100];
  int i_Cvir;
  char CvirFile[100];
} Files;

/* Various tasks the the program will perform
 */
extern struct perform_tasks {
  int All;
  int real_space_xi;
  int z_space_xi;
  int kaiser_xi;
  int multipoles;
  int r_half;
  int wp_minimize;
  int zspace_minimize;
  int MCMC;
  int HOD;
  int PVD;
  int populate_sim;
  int matter_xi;
  int matter_pk;
  int sigma_r;
  int cvir;
  char root_filename[100];
} Task;


