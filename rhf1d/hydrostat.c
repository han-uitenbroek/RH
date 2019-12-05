/* ------- file: -------------------------- hydrostat.c -------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Feb 24 15:30:40 2011 --

       --------------------------                      ----------RH-- */

/* --- Solve for hydrostatic equilibrium --            -------------- */

 
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "background.h"
#include "constant.h"
#include "accelerate.h"
#include "error.h"


/* --- Acceleration parameters --                      -------------- */

#define NG_HSE_DELAY   0
#define NG_HSE_ORDER   0
#define NG_HSE_PERIOD  0


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry geometry;
extern char   messageStr[];


/* ------- begin -------------------------- Hydrostatic.c ----------- */

void Hydrostatic(int NmaxIter, double iterLimit)
{
  const char routineName[] = "Hydrostatic";
  register int k, n, m, i;

  bool_t  Debeye, quiet, H2present;
  int    niter, Ngdelay, Ngperiod, Ngorder, Nhse;
  double *n_k, **dfdn, C1, C2, C3, beta, Phi_H, Phi_H2, *F, *dFdne, *f,
    dnmax, nHtot_old, *np, *nH2;
  struct  Ng *Nghse;
  Atom *atom;

  C1 = AMU / (2.0 * KBOLTZMANN);
  C2 = atmos.gravity / KBOLTZMANN;
  if (atmos.Stokes) C3 = 1.0 / (2.0 * MU_0 * KBOLTZMANN);

  F     = (double *) malloc(atmos.Nspace * sizeof(double));
  dFdne = (double *) malloc(atmos.Nspace * sizeof(double));
  
  /* --- Get the electron fraction F and its numerical derivative.
         ne^{met} = F(T, ne) n_H^{tot} --              -------------- */

  FMetals(F);
  dFMetals(dFdne);

  /* --- Go through all depth points and iterate on the number
         densities of nH, ne, np, nH2, and nHtot --    -------------- */

  np  = atmos.H->n[atmos.H->Nlevel-1];
  nH2 = atmos.H2->n;

  for (k = 0; k < atmos.Nspace;  k++) {

    if (atmos.T[k] >= atmos.molecules[0].Tmin &&
	atmos.T[k] <= atmos.molecules[0].Tmax) {
      H2present = TRUE;
      Nhse = 5;
    } else {
      H2present = FALSE;
      Nhse = 4;
    }
    f    = (double *) malloc(Nhse * sizeof(double));
    n_k  = (double *) malloc(Nhse * sizeof(double)); 
    dfdn = matrix_double(Nhse, Nhse);

    /* --- Get the equilibrium constants for hydrogen ionization and
           H2 association and dissociation. This can be done because we
           keep temperature and radiation field constant during the
           hydrostatic equilibrium iterations --       -------------- */

    beta = (atmos.totalAbund - 1.0) + 
      atmos.wght_per_H * C1 * SQ(atmos.vturb[k]) / atmos.T[k];

    if (H2present)
      Phi_H2 = nH2[k] / SQ(atmos.H->ntotal[k]);

    /* --- Starting solution:
           Note: n_k[0] is the amount of atomic hydrogen. -- -------- */

    n_k[0] = atmos.H->ntotal[k] - np[k];
    n_k[1] = atmos.ne[k];
    n_k[2] = np[k];
    if (H2present) n_k[3] = nH2[k];
    n_k[Nhse-1] = atmos.nHtot[k];

    Phi_H  = n_k[0] / (n_k[1] * n_k[2]);

    /* --- Initialize structure for Ng convergence acceleration -- -- */

    Nghse = NgInit(Nhse, Ngdelay=NG_HSE_DELAY,
		   Ngorder=NG_HSE_ORDER, Ngperiod=NG_HSE_PERIOD, n_k);

    niter = 1;
    while (niter <= NmaxIter) {

      for (n = 0;  n < Nhse;  n++) {
	f[n] = 0.0;
	for (m = 0;  m < Nhse;  m++) dfdn[n][m] = 0.0;
      }
      /* --- Row for mass conservation --              -------------- */

      for (m = 0;  m < Nhse-1;  m++) {
	f[0] += n_k[m];
	dfdn[0][m] = 1.0;
      }
      f[0] += beta * n_k[Nhse-1];
      dfdn[0][Nhse-1] = beta;

      /* --- Row for charge conservation --            -------------- */

      f[1] = n_k[1] - n_k[2] - F[k] * n_k[Nhse-1];
      dfdn[1][1] = 1.0 - n_k[Nhse-1] * dFdne[k];
      dfdn[1][2] = -1.0;
      dfdn[1][Nhse-1] = -F[k];

      /* --- Row for hydrogen number conservation --   -------------- */

      f[2] = n_k[0] + n_k[2] - n_k[Nhse-1];
      dfdn[2][2] = dfdn[2][0] = 1.0;
      dfdn[2][Nhse-1] = -1.0;

      if (H2present) {
        f[2] += n_k[3];
        dfdn[2][3] = 1.0;

	/* --- Row for H2 equilibrium --               -------------- */

	f[3] = -SQ(n_k[0]) * Phi_H2 + n_k[3];
	dfdn[3][0] = -2.0*n_k[0] * Phi_H2;
	dfdn[3][3] = 1.0;
      }

      /* --- Row for hydrogen ionization --            -------------- */

      f[Nhse-1] = n_k[0] - n_k[1] * n_k[2] * Phi_H;
      dfdn[Nhse-1][0] = 1.0;
      dfdn[Nhse-1][1] = -n_k[2] * Phi_H;
      dfdn[Nhse-1][2] = -n_k[1] * Phi_H;

      /* --- Fill right hand side --                   -------------- */

      f[0] -= C2 * geometry.cmass[k] / atmos.T[k];
      if (atmos.Stokes) f[0] += C3 * SQ(atmos.B[k]) / atmos.T[k];

      /* --- Solve linearized set --                   -------------- */

      SolveLinearEq(Nhse, dfdn, f, TRUE);
      for (n = 0;  n < Nhse;  n++)  n_k[n] -= f[n];

      /* --- Check convergence and accelerate if appropriate -- ----- */

      Accelerate(Nghse, n_k);
      sprintf(messageStr,
              "\n%s-- Hydrostatic equilibrium: depth %3.3d, iteration %d",
              (niter == 1) ? "\n" : "", k, niter);

      if ((dnmax = MaxChange(Nghse, messageStr, quiet=TRUE)) <= iterLimit)
        break;
      niter++;
    }
    if (dnmax > iterLimit) {
      sprintf(messageStr,
	      "Hydrostatic equilibrium iteration not converged:\n"
              " temperature: %6.1f [K], \n"
              " hydrogen density: %9.3E [m^-3],\n dnmax: %9.3E\n",
              atmos.T[k], atmos.nHtot[k], dnmax);
      Error(WARNING, routineName, messageStr);
    }
    /* --- Store results --                            -------------- */

    nHtot_old = atmos.H->ntotal[k];

    atmos.H->ntotal[k] = n_k[0] + n_k[2];
    atmos.ne[k]        = n_k[1];
    np[k]              = n_k[2];
    atmos.nHtot[k]     = n_k[Nhse-1];

    if (H2present)
      nH2[k] = n_k[3];
    else
      nH2[k] = 0.0;

    /* --- Adjust Non-LTE hydrogen level popultions -- -------------- */

    for (i = 0;  i < atmos.H->Nlevel;  i++)
      atmos.H->n[i][k] *= atmos.H->ntotal[k] / nHtot_old;   

    NgFree(Nghse);
    free(f);
    free(n_k);
    freeMatrix((void**) dfdn);    
  }
  /* --- Adjust total populations of background metals and recalculate
         all the LTE populations, including hydrogen -- ------------- */

  LTEpops(atmos.H, Debeye=TRUE);

  for (n = 1;  n < atmos.Natom;  n++) {
    atom = &atmos.atoms[n];
    for (k = 0;  k < atmos.Nspace;  k++)
      atom->ntotal[k] = atom->abundance * atmos.nHtot[k];

    LTEpops(atom, Debeye=TRUE);
  }

  /* --- Clean up --                                   -------------- */

  free(F);
  free(dFdne);
}
/* ------- end ---------------------------- Hydrostatic.c ----------- */
