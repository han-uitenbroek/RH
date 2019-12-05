/* ------- file: -------------------------- ltepops.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Nov 17 10:07:34 2010 --

       --------------------------                      ----------RH-- */

/* --- Various routines to calculate LTE populations -- ------------- */
 
#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "error.h"
#include "statistics.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern char messageStr[];


/* ------- begin -------------------------- LTEpops.c --------------- */

void LTEpops(Atom *atom, bool_t Debeye)
{
  register int k, i, m;

  char    labelStr[MAX_LINE_SIZE];
  int     Z, dZ, *nDebeye, Nspace = atmos.Nspace;
  double  cNe_T, dE_kT, dEion, dE, gi0, c1, sum, c2;

  /* --- Computes LTE populations of a given atom.
         Takes account of Debeye shielding and lowering of the ionization
         potential when Debeye is set to TRUE:

           dE_ion = -Z * (e^2/4PI*EPSILON_0) / D           [J]
                D = sqrt(EPSILON_0/(2e^2)) * sqrt(kT/ne)   [m]
 
    See: Mihalas (78), pp. 293-295

       --                                              -------------- */

  getCPU(3, TIME_START, NULL);

  /* --- Depth-independent constants --               --------------- */

  c1 = (HPLANCK/(2.0*PI*M_ELECTRON)) * (HPLANCK/KBOLTZMANN);

  /* --- Determine the total lowering of ionization potential due
         to Debeye shielding --                       --------------- */

  if (Debeye) {
    c2 = sqrt(8.0*PI/KBOLTZMANN) *
      pow(SQ(Q_ELECTRON)/(4.0*PI*EPSILON_0), 1.5);
    nDebeye = (int *) malloc(atom->Nlevel * sizeof(int));
    for (i = 1;  i < atom->Nlevel;  i++) {
      nDebeye[i] = 0;
      Z = atom->stage[i];
      for (m = 1;  m <= (atom->stage[i] - atom->stage[0]);  m++, Z++)
	nDebeye[i] += Z;
    }
  }
  /* --- Solve Saha-Boltzmann equilibrium equations --  ------------- */

  for (k = 0;  k < Nspace;  k++) {
    if (Debeye) dEion = c2 * sqrt(atmos.ne[k] / atmos.T[k]);
    cNe_T = 0.5*atmos.ne[k] * pow(c1/atmos.T[k], 1.5);
    sum   = 1.0;

    for (i = 1;  i < atom->Nlevel;  i++) {
      dE  = atom->E[i] - atom->E[0];
      gi0 = atom->g[i] / atom->g[0];
      dZ  = atom->stage[i] - atom->stage[0];

      if (Debeye)
	dE_kT = (dE - nDebeye[i] * dEion) / (KBOLTZMANN * atmos.T[k]);
      else
	dE_kT = dE / (KBOLTZMANN * atmos.T[k]);
      
      atom->nstar[i][k] = gi0 * exp(-dE_kT);
      for (m = 1;  m <= dZ;  m++) atom->nstar[i][k] /= cNe_T;
      sum += atom->nstar[i][k];
    }
    atom->nstar[0][k] = atom->ntotal[k] / sum;

    for (i = 1;  i < atom->Nlevel;  i++)
      atom->nstar[i][k] *= atom->nstar[0][k];
  }

  if (Debeye) free(nDebeye);

  sprintf(labelStr, "LTEpops %2s", atom->ID);
  getCPU(3, TIME_POLL, labelStr);
}
/* ------- end ---------------------------- LTEpops.c --------------- */

/* ------- begin -------------------------- LTEpops_elem.c ---------- */

void LTEpops_elem(Element *element)
{
  register int k, i;

  bool_t  hunt;
  double *Uk, *Ukp1, C1, *sum, *CT_ne;

  getCPU(4, TIME_START, NULL);

  C1 = (HPLANCK/(2.0*PI*M_ELECTRON)) * (HPLANCK/KBOLTZMANN);

  sum   = (double *) malloc(atmos.Nspace * sizeof(double));
  CT_ne = (double *) malloc(atmos.Nspace * sizeof(double));
  Uk    = (double *) malloc(atmos.Nspace * sizeof(double));
  Ukp1  = (double *) malloc(atmos.Nspace * sizeof(double));

  for (k = 0;  k < atmos.Nspace;  k++) {
    CT_ne[k] = 2.0 * pow(C1/atmos.T[k], -1.5) / atmos.ne[k];
    sum[k]   = 1.0;
    element->n[0][k] = 1.0;
  }

  Linear(atmos.Npf, atmos.Tpf, element->pf[0],
         atmos.Nspace, atmos.T, Uk, hunt=TRUE);

  for (i = 1;  i < element->Nstage;  i++) {
    Linear(atmos.Npf, atmos.Tpf, element->pf[i],
           atmos.Nspace, atmos.T, Ukp1, hunt=TRUE);

    for (k = 0;  k < atmos.Nspace;  k++) {
      element->n[i][k] = element->n[i-1][k] * CT_ne[k] *
        exp(Ukp1[k] - Uk[k] -
	    element->ionpot[i-1]/(KBOLTZMANN*atmos.T[k]));
      sum[k] += element->n[i][k];
    }
    SWAPPOINTER(Uk, Ukp1);
  }

  for (k = 0;  k < atmos.Nspace;  k++)
    element->n[0][k] = element->abund * atmos.nHtot[k] / sum[k];
  for (i = 1;  i < element->Nstage;  i++) {
    for (k = 0;  k < atmos.Nspace;  k++)
      element->n[i][k] *= element->n[0][k];
  }

  free(sum);     free(CT_ne);
  free(Uk);      free(Ukp1);
}
/* ------- end ---------------------------- LTEpops_elem.c ---------- */

/* ------- begin -------------------------- LTEmolecule.c ----------- */

void LTEmolecule(Molecule *molecule)
{
  /* --- Calculate partition functions for each molecular vibrational
         state v of the molecule. LTE populations are then given by:

         nv^*[k] = molecule->n * pfv[v][k] / pf[k].

   Note: The actual LTE populations are calculated (in initSolution)
         only after chemical equilibrium has been established.
         --                                            -------------- */

  register int k, v, J, kr;

  char    labelStr[MAX_LINE_SIZE];
  double  kT, gJ, **E;
  MolecularLine *mrt;

  if (!molecule->active) {
    sprintf(messageStr, "Molecule must be active: %s\n", molecule->ID);
    Error(ERROR_LEVEL_2, "LTEmolecule", messageStr);
  }
  /* --- Fill energy matrix --                         -------------- */

  E = matrix_double(molecule->Nv, molecule->NJ);
  for (kr = 0;  kr < molecule->Nrt;  kr++) {
    mrt = molecule->mrt + kr;

    E[mrt->vi][(int) (mrt->gi - 1)/2] = mrt->Ei;
    E[mrt->vj][(int) (mrt->gj - 1)/2] = mrt->Ej;
  }

  for (k = 0;  k < atmos.Nspace;  k++)
    molecule->pf[k] = 0.0;

  for (v = 0;  v < molecule->Nv;  v++) {
    for (J = 0;  J < molecule->NJ;  J++) {
      gJ = 2*J + 1;
      for (k = 0;  k < atmos.Nspace;  k++)
	molecule->pfv[v][k] +=
	  gJ * exp(-E[v][J] / (KBOLTZMANN * atmos.T[k]));
    }
    /*  --- Also store the total partition function here -- --------- */

    for (k = 0;  k < atmos.Nspace;  k++)
      molecule->pf[k] += molecule->pfv[v][k];
  }

  freeMatrix((void **) E);

  sprintf(labelStr, "LTEpops %3s", molecule->ID);
  getCPU(4, TIME_POLL, labelStr);
}
/* ------- end ---------------------------- LTEmolecule.c ----------- */

/* ------- begin -------------------------- SetLTEQuantities.c ------ */

void SetLTEQuantities(void)
{
  register int n;

  bool_t Debeye = TRUE;
  Atom *atom;

  for (n = 0;  n < atmos.Natom;  n++) {
    atom = &atmos.atoms[n];

    /* --- Get LTE populations for each atom --        -------------- */

    LTEpops(atom, Debeye);

    if (atom->active) {
      
      /* --- Read the collisional data (in MULTI's GENCOL format).
             After this we can close the input file for the active
             atom. --                                  -------------- */

      CollisionRate(atom, atom->fp_input);

      /* --- Compute the fixed rates and store in Cij -- ------------ */

      if (atom->Nfixed > 0) FixedRate(atom);
    }
  }
}
/* ------- end ---------------------------- SetLTEQuantities.c ------ */

