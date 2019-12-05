/* ---------------------------------------- Iterate.c --------------- */
/* ------- One-D, plane-parallel version ---------------------------- */

/* Implements a simple single-line scattering problem with given
 * photon destruction probability epsilon and Planck function Bp.
 *
 * Han Uitenbroek
 * Last modified: Thu May 31 09:19:23 2018 --
 */
 
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "accelerate.h"
#include "inputs.h"
#include "error.h"
#include "statistics.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern bool_t shortchar;
extern int Nlambda;
extern double *Bp, *epsilon, *phi, *wlamb, wphi, *chi, *Sny, **Iemerge;
extern char    messageStr[];
extern struct  Ng *NgS;

extern Geometry geometry;
extern InputData input;
extern char    messageStr[];


/* ------- begin -------------------------- Iterate.c --------------- */

void Iterate(int NmaxIter, double iterLimit)
{
  const char routineName[] = "Iterate";
  register int mu, la, l, k, nIter, to_obs;

  int    Ndep = geometry.Ndep, Nrays = geometry.Nrays, Nspace, lamu;
  bool_t P_only, quiet, accel;
  enum   FeautrierOrder F_order;
  double wqmu, *dJny, *chi_la, *diagonal, *P, **I, *Psi, **Psi_sc;

  FILE *fp_iter;

  dJny     = (double *) malloc(Ndep * sizeof(double));
  chi_la   = (double *) malloc(Ndep * sizeof(double));
  diagonal = (double *) malloc(Ndep * sizeof(double));

  Iemerge  = matrix_double(Nrays, Nlambda);

  if (shortchar) {
    Nspace = Ndep;
    I      = matrix_double(2*Nrays, Nspace);
    Psi_sc = matrix_double(2*Nrays, Nspace);
  } else {
    P   = (double *) malloc(Ndep * sizeof(double));
    Psi = (double *) malloc(Ndep * sizeof(double));
  }
  /* --- Write iterations  to file --                 --------------- */

  fp_iter = fopen("iter_1d_Snu", "w");

  /* --- Start main iteration loop --                 --------------- */

  nIter = 1;
  while (nIter <= NmaxIter) {
    getCPU(2, TIME_START, NULL);

    for (l = 0;  l < Ndep;  l++) {
      dJny[l] = diagonal[l] = 0.0;
    }
    if (shortchar) {
      for (la = 0;  la < Nlambda;  la++) {
	for (k = 0;  k < Nspace;  k++)
	  chi_la[k] = phi[la] * chi[k];

	/* --- Solve transfer equation for all angles -- -------------- */

	for (mu = 0;  mu < Nrays;  mu++) {
	  for (to_obs = 0;  to_obs <= 1;  to_obs++) {
	    lamu = 2*mu + to_obs;

	    if (input.S_interpolation == S_LINEAR) {
	      Piecewise_Linear_1D(la, mu, to_obs, chi_la, Sny,
				  I[lamu], Psi_sc[lamu]);
	    } else if (input.S_interpolation == S_PARABOLIC) {
	      Piecewise_1D(la, mu, to_obs, chi_la, Sny, I[lamu],
			   Psi_sc[lamu]);
	    } else if (input.S_interpolation == S_BEZIER3) {
	      Piecewise_Bezier3_1D(la, mu, to_obs, chi_la, Sny,
				   I[lamu], Psi_sc[lamu]);
	    } else {
	      sprintf(messageStr,
		      "Unknown radiation solver: %d",
		      input.S_interpolation);
	      Error(ERROR_LEVEL_1, routineName, messageStr);
	    }
	    if (to_obs) Iemerge[mu][la] = I[lamu][0];
	  }
	}
	/* --- Fill lambda operator --                   -------------- */

	for (mu = 0, lamu = mu;  mu < Nrays;  mu++) {
	  wqmu = 0.5*geometry.wmu[mu] * phi[la]*wlamb[la];
	  for (to_obs = 0;  to_obs <= 1;  to_obs++, lamu++) {
	    for (k = 0;  k < Nspace;  k++) {
	      dJny[k]     += wqmu * I[lamu][k];
	      diagonal[k] += wqmu * Psi_sc[lamu][k];
	    }
	  }
	}
      }
    } else {
      for (mu = 0;  mu < Nrays;  mu++) {
	for (la = 0;  la < Nlambda;  la++) {
	  for (l = 0;  l < Ndep;  l++) {
	    chi_la[l] = phi[la] * chi[l];
	  }

	  /* --- Formal solution and diagonal operator - ------------ */

	  Iemerge[mu][la] =
	    Feautrier(la, mu, chi_la, Sny, F_order=STANDARD, P, Psi);

	  /* --- Angle and wavelength integration --    ------------- */

	  wqmu = geometry.wmu[mu] * phi[la]*wlamb[la];
	  for (l = 0;  l < Ndep;  l++) {
	    dJny[l]     += wqmu*P[l];
	    diagonal[l] += wqmu*Psi[l];
	  }
	}
      }
    }
    /* --- Accelerated lambda iteration --              ------------- */

    for (l = 0;  l < Ndep;  l++) {
      diagonal[l] *= wphi;
      dJny[l] = wphi*dJny[l] - diagonal[l] * Sny[l];
      Sny[l]  = ((1.0 - epsilon[l])*dJny[l] + epsilon[l]*Bp[l]) /
        (1.0 - (1.0 - epsilon[l])*diagonal[l]);
    }
    /* --- Ng acceleration if requested --              ------------- */

    accel = Accelerate(NgS, Sny);
    sprintf(messageStr, "-- Main iteration No: %d", nIter);
    if (MaxChange(NgS, messageStr, quiet=FALSE) <= iterLimit)  break;
    Error(MESSAGE, NULL, (accel == TRUE) ? " (accelerated)\n" : "\n");

    fwrite(Sny, sizeof(double), Ndep, fp_iter);

    getCPU(2, TIME_POLL, "Iterate");
    nIter++;
  }

  free(dJny);  free(chi_la);  free(diagonal);
  if (shortchar) {
    freeMatrix((void **) I);
    freeMatrix((void **) Psi_sc);
  } else {
    free(P);
    free(Psi);
  }
  NgFree(NgS);

  fclose(fp_iter);
}
/* ------- end ---------------------------- Iterate.c --------------- */
