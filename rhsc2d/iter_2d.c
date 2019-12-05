/* ------- file: -------------------------- iter_2d.c --------------- */
/* ------- Two-D, plane-parallel version ---------------------------- */

/* Implements a simple single-frequency scattering problem with given
 * photon destruction probability epsilon and Planck function Bp.
 *
 * Han Uitenbroek
 * Last modified: Thu May 31 11:20:46 2018 --
 */
 
#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "inputs.h"
#include "error.h"
#include "accelerate.h"
#include "statistics.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry geometry;
extern InputData input;

extern    int  Nlambda;
extern double *Bp, *chi, *S, **Iemerge, *epsilon, *phi, *wlamb, wphi;
extern struct  Ng *Ng_S;
extern char    messageStr[];

/* ------- begin -------------------------- Iterate.c --------------- */

void Iterate(int NmaxIter, double iterLimit)
{
  const char routineName[] = "Iterate";
  register int mu, la, k, to_obs;

  bool_t  quiet, accel;
  int     niter, Nspace = atmos.Nspace, Nrays = atmos.Nrays, lamu;
  double *dJny, *diagonal, **I, **Psi, *chi_la, dSmax, wqmu;

  getCPU(1, TIME_START, NULL);

  I        = matrix_double(2*Nrays, Nspace);
  Psi      = matrix_double(2*Nrays, Nspace);

  dJny     = (double *) malloc(Nspace * sizeof(double));
  diagonal = (double *) malloc(Nspace * sizeof(double));
  chi_la   = (double *) malloc(Nspace * sizeof(double));

  Iemerge  = matrix_double(Nrays, Nlambda);

  niter = 1;
  while (niter <= NmaxIter) {
    getCPU(2, TIME_START, NULL);
    for (k = 0;  k < Nspace;  k++) {
      dJny[k] = diagonal[k] = 0.0;
    }
    for (la = 0;  la < Nlambda;  la++) {
      for (k = 0;  k < Nspace;  k++) chi_la[k] = phi[la] * chi[k];

      /* --- Solve transfer equation for all angles -- -------------- */

      for (mu = 0;  mu < Nrays;  mu++) {
	for (to_obs = 0;  to_obs <= 1;  to_obs++) {
	  lamu = 2*mu + to_obs;

	  if (input.S_interpolation == S_LINEAR) {
	    Piecewise_Linear_2D(&geometry, la, mu, to_obs,
				chi_la, S, I[lamu], Psi[lamu]);
	  } else if (input.S_interpolation == S_PARABOLIC) {
	    Piecewise_2D(&geometry, la, mu, to_obs,
			 chi_la, S, I[lamu], Psi[lamu]);
	  } else if (input.S_interpolation == S_BEZIER3) {
	    Piecewise_Bezier3_2D(&geometry, la, mu, to_obs,
				chi_la, S, I[lamu], Psi[lamu]);
	  } else {
	    sprintf(messageStr,
		    "Unknown radiation solver: %d",
		    input.S_interpolation);
	    Error(ERROR_LEVEL_1, routineName, messageStr);
	  }
	}
      }
      /* --- Fill lambda operator --                   -------------- */

      for (mu = 0, lamu = mu;  mu < Nrays;  mu++) {
        wqmu = 0.5*geometry.wmu[mu] * phi[la]*wlamb[la];
	for (to_obs = 0;  to_obs <= 1;  to_obs++, lamu++) {
	  for (k = 0;  k < Nspace;  k++) {
	    dJny[k]     += wqmu * I[lamu][k];
	    diagonal[k] += wqmu * Psi[lamu][k];
	  }
	}
      }
    }
    /* --- Accelerated lambda iteration --              ------------- */
 
    for (k = 0;  k < Nspace;  k++) {
      diagonal[k] *= wphi;
      dJny[k] =  wphi*dJny[k] - diagonal[k] * S[k];
      S[k] = ((1.0 - epsilon[k])*dJny[k] + epsilon[k]*Bp[k]) /
	(1.0 - (1.0 - epsilon[k])*diagonal[k]);
    }
    accel = Accelerate(Ng_S, S);
    sprintf(messageStr, "-- Main iteration: %3d,", niter);
    dSmax = MaxChange(Ng_S, messageStr, quiet=FALSE);
    Error(MESSAGE, NULL, (accel) ? " (accelerated)\n" : "\n");

    sprintf(messageStr, "Total Iteration %3d", niter);
    getCPU(2, TIME_POLL, messageStr);

    if (dSmax < iterLimit) break;
    niter++;
  }

  freeMatrix((void **) I);  
  freeMatrix((void **) Psi);

  free(dJny);  free(diagonal);  free(chi_la);
  NgFree(Ng_S);

  getCPU(1, TIME_POLL, "Iteration Total");
}
/* ------- end ---------------------------- Iterate.c --------------- */
