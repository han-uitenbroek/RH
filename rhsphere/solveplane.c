/* ------- file: -------------------------- solveplane.c ------------

       Version:       rh2.0, 1-D spherically symmetric
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Oct 16 11:03:34 2000   --

       --------------------------                      ----------RH-- */

/* --- Obtain 1-D plane solution with Eddington approximation (one angle)
       and Feautrier formal solution --                -------------- */

 
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "accelerate.h"
#include "constant.h"
#include "inputs.h"
#include "error.h"
#include "statistics.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- Iterate.c --------------- */

void solve_1D_plane(struct Atom *atom, Molecule *molecule,
		    int NmaxIter, double iterLimit)
{
  const char routineName[] = "solve_1D_plane";
  register int k, n, nspect, niter;
  
  bool_t   quiet, includesBB, to_obs, initialize, crosscoupling,
           boundbound, polarized, PRD_angle_dep, solve_Stokes, accel;
  enum     FeautrierOrder F_order;
  int      Nspace = atmos.Nspace, mu;
  double   dPopsmax, *chi, *P, *Psi, *J, *Jdag, *S, xmu, wmu,
          *chi_as, *eta_as, vboundary[4] = {0.0, 0.0, 0.0, 0.0};
  ActiveSet *as;

  /* --- There is only one direction in this case --   -------------- */

  mu  = 0;
  xmu = 1.0 / sqrt(3.0);
  wmu = 1.0;

  /* --- Allocate space for temporary variables like the up- and
         downward intensity, and associated diagonal operator Psi
         the opacity, emissivity and source function, and finally the
         the rate coefficient matrix Gamma. --            ----------- */

  chi    = (double *) malloc(Nspace * sizeof(double));
  S      = (double *) malloc(Nspace * sizeof(double));
  P      = (double *) malloc(Nspace * sizeof(double));
  Psi    = (double *) malloc(Nspace * sizeof(double));
  J      = (double *) malloc(Nspace * sizeof(double));
  Jdag   = (double *) malloc(Nspace * sizeof(double));
  chi_as = (double *) malloc(Nspace * sizeof(double));
  eta_as = (double *) malloc(Nspace * sizeof(double));

  /* --- Start of the main iteration loop --             ------------ */

  niter = 1;
  while (niter <= NmaxIter) {
 
    /* --- zero the radiative rates and add the collisional and
           fixed transition rates --                     ------------ */
 
    zeroRates(atom);
    initGamma(atom);
    initGammaMolecule(molecule);
 
    /* --- Loop over the spectrum to obtain formal solutions,
           get radiative contributions to coefficient matrix Gamma,
           and to invert the approximate operator --     ------------ */

    getCPU(3, TIME_START, NULL);
    for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {

      /* --- Retrieve active set as of transitions and read current
             mean intensities at wavelength nspect. --   ------------ */

      alloc_as(nspect, crosscoupling=TRUE);
      as = &spectrum.as[nspect];
      for (k = 0;  k < Nspace;  k++)  J[k] = 0.0;
      readJlambda(nspect, Jdag);

      /* --- Check whether current active set includes a bound-bound
             and/or polarized transition and/or angledependent PRD
             transition. Otherwise, only angle-independent opacity and
             source functions are needed --            -------------- */ 

      boundbound    = containsBoundBound(as);
      PRD_angle_dep = (containsPRDline(as) && input.PRD_angle_dep);
      polarized     = containsPolarized(as);
      solve_Stokes  = (polarized && input.StokesMode == FULL_STOKES);

      /* --- Case of angle-dependent opacity and source function - -- */

      if (polarized || PRD_angle_dep ||
	  (atmos.moving && (boundbound || atmos.backgrBB[nspect]))) {
	Error(ERROR_LEVEL_2, routineName,
	      "Angle dependent source functions and opacities\n"
	      " not yet implemented in spherical symmetry");
      } else {

        /* --- Angle-independent case --               -------------- */

	Opacity(nspect, 0, 0, chi_as, eta_as, initialize=TRUE);
	addtoCoupling(nspect);
	readBackground(nspect, 0, 0);

	for (k = 0;  k < Nspace;  k++) {
	  chi[k] = chi_as[k] + atmos.chi_c[k];
	  S[k]   = (eta_as[k] + atmos.eta_c[k] +
		    atmos.scatt[k]*Jdag[k]) / chi[k];
	}
	vboundary[1] = geometry.Itop[nspect];
	vboundary[3] = geometry.Icore[nspect];
	spectrum.I[nspect][mu] =
	  Feautrier(Nspace, geometry.r, xmu,
		    chi, S, vboundary, F_order=STANDARD, P, TRUE);
	Diagonal(Nspace, Psi, F_order=STANDARD);

	for (k = 0;  k < Nspace;  k++) {
	  Psi[k] /= chi[k];
	  J[k]   += P[k] * wmu;
	}
	addtoRates(nspect, wmu, P);
	addtoGamma(nspect, wmu, eta_as, P, Psi);
      }
      /* --- Write new J for current wavelength --     -------------- */

      writeJlambda(spectrum.fp_J, nspect, J);
      free_as(nspect, crosscoupling=TRUE);
    }

    if (atom) {
      statEquil(atom, input.isum);
 
      accel = Accelerate(atom->Ng_n, atom->n[0]);
      sprintf(messageStr, "\niteration #%d,", niter);
      dPopsmax = MaxChange(atom->Ng_n, messageStr, quiet=FALSE);
      Error(MESSAGE, NULL, (accel == TRUE) ? " (accelerated)\n" : "\n");
 
      /* --- Redistribute intensity in PRD lines if necessary -- ---- */
 
      if ((atom->Nprd > 0) && (dPopsmax < input.PRD_treshold))
        Redistribute(atom, input.PRD_NmaxIter, MAX(dPopsmax, iterLimit));
    } else if (molecule) {
      statEquilMolecule(molecule, 0);
 
      accel = Accelerate(molecule->Ng_nv, molecule->nv[0]);
      sprintf(messageStr, "\niteration #%d,", niter);
      dPopsmax = MaxChange(molecule->Ng_nv, messageStr, quiet=FALSE);
      Error(MESSAGE, NULL, (accel == TRUE) ? " (accelerated)\n" : "\n");
    }

    getCPU(3, TIME_POLL, "Solve Plane");
    if (dPopsmax < iterLimit) break;
    niter++;
  }
  /* --- Free space for temporary variables --          ------------- */

  free(P);       free(Psi);      free(J);       free(Jdag);
  free(S);       free(chi);      free(chi_as);  free(eta_as);
}
/* ------- end ---------------------------- solve_15D.c ------------- */
