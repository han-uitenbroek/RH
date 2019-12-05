/* ------- file: -------------------------- metal.c -----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Jul  7 16:01:45 2009 --

       --------------------------                      ----------RH-- */

/* --- Bound-free and bound-bound opacity and emissivity due to elements
       in the array metals (struct Atom *). Also opacity due to hydrogen
       bound-bound transitions is calculated in passive_bb.

       The routines duplicateLevel(active_atom, labeli) and
       duplicateLine(active_atom, labeli, labelj) are used to check
       whether the level labeli, or radiative transition between labeli
       and labelj, respectively, are part of the active set of
       transitions. In this way we try to prevent counting opacity twice,
       both in the active and the background transition.

       Global variables:
          atmos -- Atmos structure for atmospheric data.
           atom -- Atom structure with the active atom.

       Input:
         lambda -- Wavelength [nm] for which opacity and emissivity 
                   are to be calculated.
         Nmetal -- Number of entries in metals array. 
         metals -- Pointer to array of Atom structures containing atomic
                   data for metals.

       Additional input for passive_bb:
         nspect -- Index of spectrum (needed to set atmos.backgrflags[]).
             mu -- Index of ray.
         to_obs -- Boolean set to TRUE is ray is followed in direction
                   towards observer.

       Output:
         chi[Nspace] -- Array for opacities [m^2].
         eta[Nspace] -- Array for emissivities [J s^-1 Hz^-1 sr^-1].

       --                                              -------------- */
 
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "constant.h"
#include "background.h"
#include "inputs.h"


#define N_MAX_OVERLAP  10


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern char messageStr[];
extern InputData input;


/* ------- begin -------------------------- Metal_bf.c -------------- */

bool_t Metal_bf(double lambda, int Nmetal, struct Atom *metals,
		double *chi, double *eta)
{  
  register int  m, k, kr;

  bool_t   hunt;
  int      i, j, Z;
  double   lambdaEdge, alpha_la, twohnu3_c2, twohc, gijk, hc_k, hc_kla,
         **n, *expla = NULL, n_eff, gbf_0;
  Atom *metal;
  AtomicContinuum *continuum;

  twohc  = (2.0 * HPLANCK * CLIGHT) / CUBE(NM_TO_M);
  hc_k   = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M);
  for (k = 0;  k < atmos.Nspace;  k++) {
    chi[k] = 0.0;
    eta[k] = 0.0;
  }

  /* --- Go through the bound-free transitions of the metals and add
         the opacity and emissivity for each transition for which the
         current wavelength falls below treshold and above the 
         minimum wavelength --                         -------------- */
  
  for (m = 0, metal = metals;  m < Nmetal;  m++, metal++) {
    if (!metal->active) {

      /* --- Use LTE or NonLTE population numbers ? -- -------------- */

      n = (metal->n != metal->nstar) ? metal->n : metal->nstar;

      for (kr = 0;  kr < metal->Ncont;  kr++) {
	continuum = metal->continuum + kr;
	i = continuum->i;
	j = continuum->j;
	lambdaEdge = continuum->lambda0;

	if (lambda <= lambdaEdge  &&  lambda >= continuum->lambda[0]) {
	  hc_kla     = hc_k / lambda;
	  twohnu3_c2 = twohc / CUBE(lambda);

	  /* --- Evaluate the exponential only once at wavelength lambda,
	     not for each transition seperately -- ------------ */

	  if (expla == NULL) {
	    expla = (double *) malloc(atmos.Nspace * sizeof(double));
	    for (k = 0;  k < atmos.Nspace;  k++)
	      expla[k] = exp(-hc_kla/atmos.T[k]);
	  }

	  if (continuum->hydrogenic) {
	    Z = metal->stage[continuum->j];
	    n_eff = Z*sqrt(E_RYDBERG / (metal->E[continuum->j] -
					metal->E[continuum->i]));
	    gbf_0 = Gaunt_bf(continuum->lambda0, n_eff, Z);

	    alpha_la = continuum->alpha0 * CUBE(lambda/continuum->lambda0) *
	      Gaunt_bf(lambda, n_eff, Z) / gbf_0;
	  } else {
	    splineCoef(continuum->Nlambda, continuum->lambda,
		       continuum->alpha);
	    splineEval(1, &lambda, &alpha_la, hunt=FALSE);
	  }

	  for (k = 0;  k < atmos.Nspace;  k++) {
	    gijk    = metal->nstar[i][k]/metal->nstar[j][k] * expla[k];
	    chi[k] += alpha_la * (1.0 - expla[k]) * n[i][k];
	    eta[k] += twohnu3_c2 * gijk * alpha_la * n[j][k];
	  }
	}
      }
    }
  }
  if (expla != NULL) {
    free(expla);
    return TRUE;
  } else
    return FALSE;
}
/* ------- end ---------------------------- Metal_bf.c -------------- */

/* ------- begin -------------------------- passive_bb.c ------------ */

flags passive_bb(double lambda, int nspect, int mu, bool_t to_obs,
		  double *chi, double *eta, double *chip)
{  
  const char routineName[] = "passive_bb";
  register int k, kr, l, m, nc;

  static bool_t initialize = TRUE;
  static int Nlist;
  static struct Linelist *linelist[N_MAX_OVERLAP];

  bool_t   add_to_list, linepresent;
  int      i, j, entry;
  double   dlambda, phi, v, twohnu3_c2, hc, fourPI, hc_4PI,
           gij, Vij, **n;
  Atom *atom;
  AtomicLine *line;
  flags backgrflags;

  /* --- Calculate contribution of bound-bound transitions in the 
         background atoms (including hydrogen) to the opacity and
         emissivity.
 
   Note: A list of lines is maintained to prevent recalculation of
         the damping parameter of the lines for successive wavelengths
         and angles.
         --                                            -------------- */

  backgrflags.hasline     = FALSE;
  backgrflags.ispolarized = FALSE;

  if (initialize) {
    for (l = 0;  l < N_MAX_OVERLAP;  l++) linelist[l] = NULL;
    Nlist = 0;
    initialize = FALSE;
  }

  hc     = HPLANCK * CLIGHT;
  fourPI = 4.0 * PI;
  hc_4PI = hc / fourPI;

  for (k = 0;  k < atmos.Nspace;  k++) {
    chi[k] = 0.0;
    eta[k] = 0.0;
  }
  if (atmos.Stokes) {
    for (k = atmos.Nspace;  k < 4*atmos.Nspace;  k++) {
      chi[k] = 0.0;
      eta[k] = 0.0;
    }
    if (input.magneto_optical)
      for (k = 0;  k < 3*atmos.Nspace;  k++) chip[k] = 0.0;
  }
  /* --- Reset the used tags in the linelist --        -------------- */

  for (l = 0;  l < Nlist;  l++) linelist[l]->used = FALSE;

  /* --- Go through the bound-bound transitions, First hydrogen, then
         the metals, and add the opacity and emissivity for each
         transition for which the current wavelength falls within the
         limits of the line. --                        -------------- */
  
  for (m = 0;  m < atmos.Natom;  m++) {
    atom = atmos.atoms + m;
    if (!atom->active) {

      /* --- Use LTE or NonLTE population numbers ? -- -------------- */

      n = (atom->n != atom->nstar) ? atom->n : atom->nstar;

      for (kr = 0;  kr < atom->Nline;  kr++) {
	line = atom->line + kr;
	i = line->i;
	j = line->j;
	dlambda = line->lambda0 *
	  line->qwing * (atmos.vmicro_char / CLIGHT);

	if (fabs(lambda - line->lambda0) <= dlambda) {
	  backgrflags.hasline = TRUE;
	  atmos.backgrflags[nspect].hasline = TRUE;

	  /* --- Add line to list if not yet present -- ----------- */ 

	  add_to_list = TRUE;
	  for (l = 0;  l < Nlist;  l++) {
	    if (line == linelist[l]->line) {
	      add_to_list = FALSE;
	      entry = l;
	      break;
	    }
	  }
	  if (add_to_list) {
	    if (Nlist == N_MAX_OVERLAP) {
	      sprintf(messageStr, "Too many overlapping transitions");
	      Error(ERROR_LEVEL_2, routineName, messageStr);
	    }
	    /* --- Create a new entry in the list -- -------------- */

	    linelist[Nlist] =
	      (struct Linelist *) malloc(sizeof(struct Linelist));
	    entry = Nlist++;
	    linelist[entry]->line = line;

	    /* --- Calculate and store the line's damping parameter */

	    if (line->Voigt) {
	      linelist[entry]->adamp =
		(double *) malloc(atmos.Nspace * sizeof(double));
	      Damping(line, linelist[entry]->adamp);
	    } else
	      linelist[entry]->adamp = NULL;
	  }
	  linelist[entry]->used = TRUE;

	  gij = line->Bji / line->Bij;
	  twohnu3_c2 = line->Aji / line->Bji;
 
	  /* --- Evaluate absorption and emission coefficients -- - */

	  for (nc = 0;  nc < line->Ncomponent;  nc++) {
	    for (k = 0;  k < atmos.Nspace;  k++) {
	      v = (lambda - line->lambda0 - line->c_shift[nc]) *
		CLIGHT / (line->lambda0 * atom->vbroad[k]);
	      if (atmos.moving) {
		if (to_obs)
		  v += vproject(k, mu) / atom->vbroad[k];
		else
		  v -= vproject(k, mu) / atom->vbroad[k];
	      }
	      if (line->Voigt)
		phi = Voigt(linelist[entry]->adamp[k], v, NULL,
			    ARMSTRONG) * line->c_fraction[nc];
	      else
		phi = exp(-SQ(v));
 
	      Vij = hc_4PI * line->Bij * phi / (SQRTPI*atom->vbroad[k]);
	      chi[k] += Vij * (n[i][k] - gij * n[j][k]);
	      eta[k] += twohnu3_c2 * gij * Vij * n[j][k];
	    }
	  }
	}
      }
    }
  }
  /* --- Remove lines from the line list that have not been used at this
         wavelength, and sort list --                  -------------- */

  for (l = 0;  l < Nlist;  l++) {
    if (!linelist[l]->used) {
      if (linelist[l]->adamp) free(linelist[l]->adamp);
      free(linelist[l]);
      linelist[l] = NULL;
    }
  }
  for (l = 0;  l < N_MAX_OVERLAP;  l++) {
    if (linelist[l] == NULL) {
      for (m = l+1;  m < N_MAX_OVERLAP;  m++) {
	if (linelist[m] != NULL) {
	  linelist[l] = linelist[m];
          linelist[m] = NULL;
	  break;
	}
      }
    }
  }
  /* --- Count number of entries in the list --        -------------- */

  Nlist = 0;
  for (l = 0;  l < N_MAX_OVERLAP;  l++) if (linelist[l]) Nlist++;

  return backgrflags;
}
/* ------- end ---------------------------- passive_bb.c ------------ */
