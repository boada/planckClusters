#include <stdlib.h>
#include "header.h"


/* These are the globale variables. When needed, each is set to the
 * default value. Those without defualt values are either calculated
 * or set in the batfile.
 */


double GAMMA=0.2,         /* Shape parameter of EBW power spectrum */
  HUBBLE=0.7,             /* Hubble in 100 km/s/Mpc (not realy used in code). */
  SIGMA_8=0.95,           /* Normalization of power spectrum */
  RHO_CRIT=2.775e11,      /* Critial mass density  in h^2 M_sun/Mpc^3 */
  SPECTRAL_INDX=1.0,      /* n_s -> P(k) = k^n_s */
  OMEGA_M=0.1,            /* Matter density */
  OMEGA_B=0.0,            /* Baryon density */
  DELTA_CRIT=1.686,       /* Critical overdensity for linear collapse */
  MSTAR,                  /* Mass scale at which sigm(M) = DELTA_CRIT */
  REDSHIFT,               /* redshift at which to calculate statistics */
  DELTA_HALO,             /* mean interior density to define dark matter halos */
  JENKINS_A=0.315,        /* Jenkins mass function--> normalization */
  JENKINS_B=0.61,         /* Jenkins mass function--> constant in exponential */
  JENKINS_C=3.8;          /* Jenkins mass function--> exponent in exponential */

int ITRANS=4,             /* Type of transfer function to be used */
  OUTPUT=0,               /* level of output to stderr */
  RESET_COSMOLOGY=0,      /* Flag to recalculate tabulated quantities (e.g. b(M), dn/dM...) */
  ERROR_FLAG=0;           /* Global notification of problem with qromo/zbrent/etc */

int ThisTask = 0;

struct file_parameters Files;
struct perform_tasks Task;
