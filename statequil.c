/* ------- file: -------------------------- statequil.c -------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Apr  1 14:09:35 2009 --

       --------------------------                      ----------RH-- */

/* --- Solves equations of statistical equilibrium with given rate
       matrix Gamma[Nlevel*Nlevel][Nspace].

       isum defines the row of the statistical equilibrium equations
       that is to be eliminated to enforce particle conservation. If
       set to -1 then at each spatial location the row with the largest
       population from the previous iteration will be eliminated -- - */

 
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "inputs.h"
#include "accelerate.h"
#include "statistics.h"
#include "error.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern InputData input;
extern char messageStr[];

/* ------- begin -------------------------- statEquil.c ------------- */

void statEquil(Atom *atom, int isum)
{
  register int i, j, ij, k;

  int    i_eliminate, Nlevel;
  double GamDiag, nmax_k, *n_k, **Gamma_k;

  getCPU(3, TIME_START, NULL);

  Nlevel = atom->Nlevel;

  /* --- Need temporary storage because Gamma has to be solved spatial
         point by spatial point while depth is normally the fastest
         running index --                              -------------- */

  n_k     = (double *) malloc(Nlevel * sizeof(double));
  Gamma_k = matrix_double(Nlevel, Nlevel);

  for (k = 0;  k < atmos.Nspace;  k++) {
    for (i = 0, ij = 0;  i < Nlevel;  i++) {
      n_k[i] = atom->n[i][k];
      for (j = 0;  j < Nlevel;  j++, ij++)
	Gamma_k[i][j] = atom->Gamma[ij][k];
    }

    if (isum == -1) {
      i_eliminate  = 0;
      nmax_k = 0.0;
      for (i = 0;  i < Nlevel;  i++) {
	if (n_k[i] > nmax_k) {
	  nmax_k = n_k[i];
	  i_eliminate = i;
	}
      }
    } else
      i_eliminate = isum;

    /* --- For each column i sum over rows to get diagonal elements - */

    for (i = 0;  i < Nlevel;  i++) {
      GamDiag = 0.0;
      Gamma_k[i][i] = 0.0;
      n_k[i] = 0.0;

      for (j = 0;  j < Nlevel;  j++) GamDiag += Gamma_k[j][i];
      Gamma_k[i][i] = -GamDiag;
    }
    /* --- Close homogeneous set with particle conservation-- ------- */

    n_k[i_eliminate] = atom->ntotal[k];
    for (j = 0;  j < Nlevel;  j++) Gamma_k[i_eliminate][j] = 1.0;

    /* --- Solve for new population numbers at location k -- -------- */

    SolveLinearEq(Nlevel, Gamma_k, n_k, TRUE);
    for (i = 0;  i < Nlevel;  i++) atom->n[i][k] = n_k[i];
  }

  free(n_k);
  freeMatrix((void **) Gamma_k);

  getCPU(3, TIME_POLL, "Stat Equil");
}
/* ------- end ---------------------------- statEquil.c ------------- */

/* ------- begin -------------------------- statEquilMolecule.c ----- */

void statEquilMolecule(struct Molecule *molecule, int isum)
{
  register int vi, vj, vij, k;

  int    i_eliminate, Nlevel;
  double GamDiag, nmax_k, *n_k, **Gamma_k;

  getCPU(3, TIME_START, NULL);

  Nlevel = molecule->Nv;

  /* --- Need temporary storage because Gamma has to be solved spatial
         point by spatial point while depth is normally the fastest
         running index --                              -------------- */

  n_k     = (double *) malloc(Nlevel * sizeof(double));
  Gamma_k = matrix_double(Nlevel, Nlevel);

  for (k = 0;  k < atmos.Nspace;  k++) {
    if (molecule->n[k] > 0.0) {
      for (vi = 0, vij = 0;  vi < Nlevel;  vi++) {
	n_k[vi] = molecule->nv[vi][k];
	for (vj = 0;  vj < Nlevel;  vj++, vij++)
	  Gamma_k[vi][vj] = molecule->Gamma[vij][k];
      }

      if (isum == -1) {
	i_eliminate  = 0;
	nmax_k = 0.0;
	for (vi = 0;  vi < Nlevel;  vi++) {
	  if (n_k[vi] > nmax_k) {
	    nmax_k = n_k[vi];
	    i_eliminate = vi;
	  }
	}
      } else
	i_eliminate = isum;

      /* --- For each column i sum over rows to get diagonal elements */

       for (vi = 0;  vi < Nlevel;  vi++) {
	GamDiag = 0.0;
	Gamma_k[vi][vi] = 0.0;
	n_k[vi] = 0.0;

	for (vj = 0;  vj < Nlevel;  vj++) GamDiag += Gamma_k[vj][vi];
	Gamma_k[vi][vi] = -GamDiag;
      }
      /* --- Close homogeneous set with particle conservation ------- */

      n_k[i_eliminate] = molecule->n[k];
      for (vj = 0;  vj < Nlevel;  vj++) Gamma_k[i_eliminate][vj] = 1.0;

      /* --- Solve for new population numbers at location k --------- */

      SolveLinearEq(Nlevel, Gamma_k, n_k, TRUE);
      for (vi = 0;  vi < Nlevel;  vi++) molecule->nv[vi][k] = n_k[vi];
    } else
      for (vi = 0;  vi < Nlevel;  vi++) molecule->nv[vi][k] = 0.0;
  }

  free(n_k);
  freeMatrix((void **) Gamma_k);

  getCPU(3, TIME_POLL, "Stat Equil");
}
/* ------- end ---------------------------- statEquilMolecule.c ----- */

/* ------- begin -------------------------- updatePopulations.c ----- */

double updatePopulations(int niter)
{
  register int nact;

  bool_t accel, quiet;
  double dpops, dpopsmax = 0.0;
  Atom *atom;
  Molecule *molecule;

  /* --- Update active atoms --                        -------------- */

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];

    statEquil(atom, input.isum);
    accel = Accelerate(atom->Ng_n, atom->n[0]);

    sprintf(messageStr, " %s,", atom->ID);
    dpops = MaxChange(atom->Ng_n, messageStr, quiet=FALSE);
    Error(MESSAGE, NULL, (accel) ? " (accelerated)\n" : "\n");

    dpopsmax = MAX(dpops, dpopsmax);
  }
  /* --- Update active molecules --                    -------------- */

  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];

    statEquilMolecule(molecule, 0);
    accel = Accelerate(molecule->Ng_nv, molecule->nv[0]);

    sprintf(messageStr, " %s ,", molecule->ID);
    dpops = MaxChange(molecule->Ng_nv, messageStr, quiet=FALSE);
    Error(MESSAGE, NULL, (accel) ? " (accelerated)\n" : "\n");

    dpopsmax = MAX(dpops, dpopsmax);
  }

  return dpopsmax;
}
/* ------- end ---------------------------- updatePopulations.c ----- */
