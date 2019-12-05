/* ------- file: -------------------------- fillgamma.c -------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Jul 24 12:34:18 2009 --

       --------------------------                      ----------RH-- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "atom.h"
#include "inputs.h"
#include "constant.h"

/* --- Routines for wavelength- and angle-integrated contributions to
       the crosscoupling, Gamma matrix, and radiative rates.

       Convention: \Gamma_ij = Gamma[i][j] represents the
                   transition j --> i, so that \Gamma_ij * n_j
		   is the rate per sec out of level j to level i.
       --                                              -------------- */


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- initGammaAtom.c --------- */

void initGammaAtom(Atom *atom)
{
  register int ij, k;

  /* --- Add the fixed rates into Gamma --             -------------- */

  for (ij = 0;  ij < SQ(atom->Nlevel);  ij++) {
    for (k = 0;  k < atmos.Nspace;  k++)
      atom->Gamma[ij][k] = atom->C[ij][k];
  }
}
/* ------- end ---------------------------- initGammaAtom.c --------- */

/* ------- begin -------------------------- initGammaMolecule.c ----- */

void initGammaMolecule(Molecule *molecule)
{
  register int ij, ji, k, vi, vj;

  /* --- Add the fixed rates into Gamma --             -------------- */

  for (vi = 0;  vi < molecule->Nv-1;  vi++) {
    vj = vi + 1;
    ij = vi*molecule->Nv + vj;
    ji = vj*molecule->Nv + vi;
    for (k = 0;  k < atmos.Nspace;  k++) {
      if (molecule->n[k]) {
	molecule->Gamma[ij][k] = molecule->C_ul[k];
	molecule->Gamma[ji][k] = molecule->C_ul[k] *
	  molecule->nvstar[vj][k] / molecule->nvstar[vi][k];
      }
    }
  }
}
/* ------- end ---------------------------- initGammaMolecule.c ----- */

/* ------- begin -------------------------- addtoGamma.c ------------ */

void addtoGamma(int nspect, double wmu, double *I, double *Psi)
{
  const char routineName[] = "addtoGamma";
  register int nact, n, k, m;

  int    i, j, ij, ji, jp, nt;
  double twohnu3_c2, twohc, wlamu, *Ieff,
        *Stokes_Q, *Stokes_U, *Stokes_V, *eta_Q, *eta_U, *eta_V;

  Atom *atom;
  AtomicLine *line;
  AtomicContinuum *continuum;
  Molecule *molecule;
  MolecularLine *mrt;
  ActiveSet *as;
  
  twohc = 2.0*HPLANCK*CLIGHT / CUBE(NM_TO_M);

  as = &spectrum.as[nspect];
  nt = nspect % input.Nthreads;

  if (containsActive(as)) {
    Ieff = (double *) malloc(atmos.Nspace * sizeof(double));

    if (input.StokesMode == FULL_STOKES  &&  containsPolarized(as)) {

      /* --- Use pointers to the bottom 3/4 of I and
             atom->rhth.eta --                         -------------- */

      Stokes_Q = I + atmos.Nspace;
      Stokes_U = I + 2*atmos.Nspace;
      Stokes_V = I + 3*atmos.Nspace;
    }
  }
  /* --- Contributions from the active transitions in atoms -- ------ */

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];

    if (as->Nactiveatomrt[nact] > 0) {
      if (input.StokesMode == FULL_STOKES  &&  containsPolarized(as)) {

	eta_Q = atom->rhth[nt].eta + atmos.Nspace;
	eta_U = atom->rhth[nt].eta + 2*atmos.Nspace;
	eta_V = atom->rhth[nt].eta + 3*atmos.Nspace;

	for (k = 0;  k < atmos.Nspace;  k++) {
	  Ieff[k] = I[k] + Stokes_Q[k] + Stokes_U[k] + Stokes_V[k] - 
	    Psi[k] * (atom->rhth[nt].eta[k] +
		      eta_Q[k] + eta_U[k] + eta_V[k]);
	}
      } else { 
	for (k = 0;  k < atmos.Nspace;  k++) {
	  Ieff[k] = I[k] - Psi[k] * atom->rhth[nt].eta[k];
	}
      }
    }

    for (n = 0;  n < as->Nactiveatomrt[nact];  n++) {
      switch (as->art[nact][n].type) {
      case ATOMIC_LINE:
	line = as->art[nact][n].ptype.line;
	i = line->i;
	j = line->j;
	twohnu3_c2 = line->Aji / line->Bji;
	break;

      case ATOMIC_CONTINUUM:
	continuum = as->art[nact][n].ptype.continuum;
	i = continuum->i;
	j = continuum->j;
	twohnu3_c2 = twohc / CUBE(spectrum.lambda[nspect]);
	break;

      default:
	sprintf(messageStr, "Invalid transition type");
	Error(ERROR_LEVEL_1, routineName, messageStr);
	twohnu3_c2 = 0.0;
      }

      if (input.Nthreads > 1) pthread_mutex_lock(&atom->Gamma_lock);
      
      ij = i*atom->Nlevel + j;
      ji = j*atom->Nlevel + i;

      for (k = 0;  k < atmos.Nspace;  k++) {
	wlamu = atom->rhth[nt].Vij[n][k] * atom->rhth[nt].wla[n][k] * wmu;

	atom->Gamma[ji][k] += Ieff[k] * wlamu;
	atom->Gamma[ij][k] += (twohnu3_c2 + Ieff[k]) *
	  atom->rhth[nt].gij[n][k] * wlamu;
      }
      /* --- Cross-coupling terms, currently only for Stokes_I -- --- */

      for (k = 0;  k < atmos.Nspace;  k++) {
	atom->Gamma[ij][k] -= atom->rhth[nt].chi_up[i][k] *
	  Psi[k]*atom->rhth[nt].Uji_down[j][k] * wmu;
      }
      /* --- If rt->i is also an upper level of another transition that
             is active at this wavelength then Gamma[ji] needs to be
             updated as well --                        -------------- */

      for (m = 0;  m < as->Nactiveatomrt[nact];  m++) {
	switch (as->art[nact][m].type) {
	case ATOMIC_LINE:     
	  jp = as->art[nact][m].ptype.line->j;
	  break;
	case ATOMIC_CONTINUUM:
	  jp = as->art[nact][m].ptype.continuum->j;
	  break;
	default:;
	}
	if (jp == i) {
	  for (k = 0;  k < atmos.Nspace;  k++) {
	    atom->Gamma[ji][k] += atom->rhth[nt].chi_down[j][k] *
	      Psi[k]*atom->rhth[nt].Uji_down[i][k] * wmu;
	  }
	}
      }
      if (input.Nthreads > 1) pthread_mutex_unlock(&atom->Gamma_lock); 
    }
  }
  /* --- Add the active molecular contributions --     -------------- */

  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];

    for (n = 0;  n < as->Nactivemolrt[nact];  n++) {
      switch (as->mrt[nact][n].type) {

      case VIBRATION_ROTATION:
	mrt = as->mrt[nact][n].ptype.vrline;
	i = mrt->vi;
	j = mrt->vj;
	twohnu3_c2 = mrt->Aji / mrt->Bji;
	break;

      default:
	sprintf(messageStr, "Invalid transition type");
	Error(ERROR_LEVEL_1, routineName, messageStr);
	twohnu3_c2 = 0.0;
      }

      if (input.Nthreads > 1) pthread_mutex_lock(&molecule->Gamma_lock);

      /* --- In case of molecular vibration-rotation transitions -- - */

      ij = i*molecule->Nv + j;
      ji = j*molecule->Nv + i;

      for (k = 0;  k < atmos.Nspace;  k++) {
	if (molecule->n[k]) {
	  wlamu = molecule->rhth[nt].Vij[n][k] *
	    molecule->rhth[nt].wla[n][k] * wmu;
	  molecule->Gamma[ji][k] += I[k] * wlamu;
	  molecule->Gamma[ij][k] += molecule->rhth[nt].gij[n][k] *
	    (twohnu3_c2 + I[k]) * wlamu;
	}
      }
      if (input.Nthreads > 1) pthread_mutex_unlock(&molecule->Gamma_lock);
    }
  }

  if (containsActive(as)) free(Ieff);
}
/* ------- end ---------------------------- addtoGamma.c ------------ */

/* ------- begin -------------------------- addtoCoupling.c --------- */

void addtoCoupling(int nspect)
{
  const  char routineName[] = "addtoCoupling";
  register int nact, n, k;

  int    i, j, nt;
  double twohnu3_c2, chicc, twohc, *n_i, *n_j;
  Atom *atom;
  AtomicLine *line;
  AtomicContinuum *continuum;
  ActiveSet *as;

  twohc = 2.0*HPLANCK*CLIGHT / CUBE(NM_TO_M);

  as = &spectrum.as[nspect];
  nt = nspect % input.Nthreads;

  /* --- Zero the cross coupling matrices --           -------------- */

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];

    for (n = 0;  n < as->Nlower[nact];  n++) {
      i = as->lower_levels[nact][n];
      for (k = 0;  k < atmos.Nspace;  k++)
	atom->rhth[nt].chi_up[i][k] = 0.0;
    }
    for (n = 0;  n < as->Nupper[nact];  n++) {
      j = as->upper_levels[nact][n];
      for (k = 0;  k < atmos.Nspace;  k++) {
	atom->rhth[nt].chi_down[j][k] = 0.0;
	atom->rhth[nt].Uji_down[j][k] = 0.0;
      }
    }
  }
  /* --- Gather terms for cross-coupling between overlapping
     transitions --                                -------------- */

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];

    for (n = 0;  n < as->Nactiveatomrt[nact];  n++) {
      switch (as->art[nact][n].type) {
      case ATOMIC_LINE:
	line = as->art[nact][n].ptype.line;
	i = line->i;
	j = line->j;
	n_i = atom->n[i];
	n_j = atom->n[j];
	twohnu3_c2 = line->Aji / line->Bji;
	break;

      case ATOMIC_CONTINUUM:
	continuum = as->art[nact][n].ptype.continuum;
	i = continuum->i;
	j = continuum->j;
	n_i = atom->n[i];
	n_j = atom->n[j];
	twohnu3_c2 = twohc / CUBE(spectrum.lambda[nspect]);
	break;

      default:
	sprintf(messageStr, "Invalid transition type");
	Error(ERROR_LEVEL_1, routineName, messageStr);
	twohnu3_c2 = 0.0;
      }
      /* --- Evaluate the cross-coupling coefficients -- ------------ */

      if (twohnu3_c2) {
	for (k = 0;  k < atmos.Nspace;  k++) {
	  chicc = atom->rhth[nt].Vij[n][k] * atom->rhth[nt].wla[n][k] *
	    (n_i[k] - atom->rhth[nt].gij[n][k]*n_j[k]);
	  atom->rhth[nt].chi_up[i][k]   += chicc;
	  atom->rhth[nt].chi_down[j][k] += chicc;

	  atom->rhth[nt].Uji_down[j][k] +=
	    twohnu3_c2 * atom->rhth[nt].gij[n][k] * atom->rhth[nt].Vij[n][k];
	}
      }
    }
  }
}
/* ------- end ---------------------------- addtoCoupling.c --------- */

/* ------- begin -------------------------- zeroRates.c ------------- */

void zeroRates(bool_t redistribute)
{
  register int kr, k, n;

  Atom *atom;

  /* --- Initialize the radiative rates for atomic transitions.

         When redistribute == TRUE only the rates of PRD lines are
         initialized.
         --                                            -------------- */

  for (n = 0;  n < atmos.Natom;  n++) {
    atom = &atmos.atoms[n];
    if (atom->active) {
      for (kr = 0;  kr < atom->Nline;  kr++) {
        if (!redistribute || (redistribute && atom->line[kr].PRD)) {
	  for (k = 0;  k < atmos.Nspace;  k++) {
	    atom->line[kr].Rij[k] = 0.0;
	    atom->line[kr].Rji[k] = 0.0;
	  }
	}
      }
      if (!redistribute) {
	for (kr = 0;  kr < atom->Ncont;  kr++) {
	  for (k = 0;  k < atmos.Nspace;  k++) {
	    atom->continuum[kr].Rij[k] = 0.0;
	    atom->continuum[kr].Rji[k] = 0.0;
	  }
	}
      }
    }
  }
}
/* ------- end ---------------------------- zeroRates.c ------------- */

/* ------- begin -------------------------- addtoRates.c ------------ */

void addtoRates(int nspect, int mu, bool_t to_obs, double wmu,
		double *I, bool_t redistribute)
{
  register int nact, n, k;

  int    la, lamu, nt;
  double twohnu3_c2, twohc, hc_4PI, Bijxhc_4PI, wlamu, *Rij, *Rji,
         up_rate, *Stokes_Q, *Stokes_U, *Stokes_V;

  ActiveSet *as;
  Atom *atom;
  AtomicLine *line;
  AtomicContinuum *continuum;
  pthread_mutex_t *rate_lock;

  /* --- Calculate the radiative rates for atomic transitions.

         When redistribute == TRUE only radiative rates of PRD lines
         are evaluated. These are needed in the iterative update of the
         emission profile ratio \rho.
         --                                            -------------- */

  twohc = 2.0*HPLANCK*CLIGHT / CUBE(NM_TO_M);

  as = &spectrum.as[nspect];
  nt = nspect % input.Nthreads;

  if (input.StokesMode == FULL_STOKES  && containsPolarized(as)){

    /* --- Use pointers to the bottom 3/4 of I and as->eta -- ------- */

    Stokes_Q = I + atmos.Nspace;
    Stokes_U = I + 2*atmos.Nspace;
    Stokes_V = I + 3*atmos.Nspace;
  }

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];

    for (n = 0;  n < as->Nactiveatomrt[nact];  n++) {
      switch (as->art[nact][n].type) {
      case ATOMIC_LINE:
	line = as->art[nact][n].ptype.line;
	if (redistribute && !line->PRD)
	  Rij = NULL;
	else {
	  Rij = line->Rij;
	  Rji = line->Rji;
	  twohnu3_c2 = line->Aji / line->Bji;

	  rate_lock = &line->rate_lock;
	}
	break;

      case ATOMIC_CONTINUUM:
	if (redistribute)
	  Rij = NULL;
	else {
	  continuum = as->art[nact][n].ptype.continuum;
	  Rij = continuum->Rij;
	  Rji = continuum->Rji;
	  twohnu3_c2 = twohc / CUBE(spectrum.lambda[nspect]);

	  rate_lock = &continuum->rate_lock;
	}
	break;
      
      default:
	Rij = NULL;
      }
      /* --- Convention: Rij is the rate for transition i -> j -- ----- */

      if (Rij != NULL) {
	if (input.Nthreads > 1) pthread_mutex_lock(rate_lock);

	for (k = 0;  k < atmos.Nspace;  k++) {
	  wlamu =
	    atom->rhth[nt].Vij[n][k] * atom->rhth[nt].wla[n][k] * wmu;
	  Rij[k] += I[k] * wlamu;
	  Rji[k] += atom->rhth[nt].gij[n][k] * (twohnu3_c2 + I[k]) * wlamu;
	}

	if (input.Nthreads > 1) pthread_mutex_unlock(rate_lock);
      }
    }
  }
}
/* ------- end ---------------------------- addtoRates.c ------------ */
