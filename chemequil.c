/* ------- file: -------------------------- chemequil.c -------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Feb  6 10:33:22 2018 --

       --------------------------                      ----------RH-- */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "background.h"
#include "accelerate.h"
#include "constant.h"
#include "error.h"
#include "statistics.h"
#include "inputs.h"


#define COMMENT_CHAR  "#"


/* --- Acceleration parameters --                      -------------- */

#define NG_CHEM_DELAY   0
#define NG_CHEM_ORDER   0
#define NG_CHEM_PERIOD  0


/* --- Evaluate chemical equilibrium for set of molecules and the complete
       and sufficient set of their constituent nuclei.

       Equations to solve: 

    1) Conservation of number density for each constituent:

        n_i + Sum_m {N_i^m * n_m} = A_i * nHtot,

       n_i are the atomic population number densities (all atoms and ions
       of one species not bound in molecules), abund_i the element abundance,
       n_m the molecular number densities of molecules containing element i,
       and N_i^m is the number of nuclei of element i in molecule m.


    2) Chemical equilibrium for each molecule (Saha):

        n_m = Prod_i {(f0_i * n_i)^N_i^m} * Phi_m(T),

        where f0_i is the fraction of atoms i in the neutral stage,
        and Phi_m(T) is the equilibrium constant for molecule m.
        If molecule n_m has a charge of +1, one of the f0_i is replaced
        by f1_i, the fraction of element i in the first ionization stage.

  Note: Molecules with other charge than 0 or +1 are currently not allowed
        and are rejected in routine readMolecule.c


    -- The Hminus population number is added in the hydrogen conservatiom
       equation, which is always the first equation.
       Hmin formation is given by:

        nHmin = ne * nH * PhiHmin,

       where PhiHmin = 1/4 * (h^2/(2PI m_e kT))^3/2 exp(Ediss/kT)


 Note: The total hydrogen number density (including protons, H-, and
       nuclei that are part of molecules like H2 and H2+ is stored in
       atmos->nHtot. The total number not in molecules or H- (i.e. the
       neutral atoms plus the protons is stored in atmos->H.ntotal.

 Note: Although the partial derivatives matrix df is mostly constant
       we have to refill it every iteration since the LU-decomposition
       used in the matrix inversion routine uses the matrix as scratch
       space. An alternative would be to copy the matrix before inversion.


 Note: The following order is used for the number densities in the
       Newton-Raphson scheme:

         [nH, n_nuclei[1..Nnuclei-1], n_molecules[0..atmos->Nmolecules]]

      The equations in the Newton-Raphson scheme are ordered as follows:

        - equations for number conservation of nucleus[0..Nnucl-1]
        - Saha equations for molecules[0..atmos->Nmolecules]

       --                                              -------------- */
 

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern InputData input; 
extern char   messageStr[];


/* ------- begin -------------------------- ChemicalEquilibrium.c --- */

void ChemicalEquilibrium(int NmaxIter, double iterLimit)
{
  const char routineName[] = "ChemicalEquilibrium";
  register int k, i, j, nu;

  char    tmpStr[13];
  bool_t  quiet;
  int     Nequation, **nucl_index, niter, Nnuclei, Ngdelay, Ngperiod,
          Ngorder, Nmaxstage;
  double *f, *a, *n, **df, *Phi, PhiHmin, fHmin, CI, dnmax = 0.0,
         *fn0, fraction, saha, *fjk, *dfjk;
  struct  Ng *Ngn;
  Atom *atom;
  Molecule *molecule;
  Element **nuclei;

  getCPU(3, TIME_START, NULL);

  /* --- Constant for Saha equation Hminus --          -------------- */

  CI = (HPLANCK/(2.0*PI*M_ELECTRON)) * (HPLANCK/KBOLTZMANN);


  /* --- Solve the chemical equilibrium equations. First,
         determine for what atoms and molecules the equilibrium
         equations have to be solved. --               -------------- */

  /* --- Collect pointers to elements that can be bound in molecules. */

  Nnuclei = 0;
  nuclei  = (Element **) malloc(atmos.Nelem * sizeof(Element *));
  for (i = 0;  i < atmos.Nelem;  i++) {
    if (atmos.elements[i].Nmolecule > 0) {
      nuclei[Nnuclei++] = &atmos.elements[i];
    }
  }
  nuclei = (Element **) realloc(nuclei, Nnuclei * sizeof(Element *));

  /* --- Check that first Nucleus is hydrogen --       -------------- */

  if (!strstr(nuclei[0]->ID, "H ")) {
    sprintf(messageStr, "First nucleus must be H not %s "
	    "(check H2.molecule)", nuclei[0]->ID);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  Nmaxstage = 0;
  for (j = 0;  j < Nnuclei;  j++) {
    if (nuclei[j]->model != NULL) {
      if (nuclei[j]->model->stage[0] > 0) {
	sprintf(messageStr,
		"Model for element %s does not have a neutral stage\n"
		" needed for molecular formation\n",
		nuclei[j]->ID);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
    } else
      Nmaxstage = MAX(Nmaxstage, nuclei[j]->Nstage);
  }
  if (Nmaxstage) {
    fjk  = (double *) malloc(Nmaxstage * sizeof(double));
    dfjk = (double *) malloc(Nmaxstage * sizeof(double));
  }

  /* --- Quantity nucl_index[i][j] stores the index (in array nuclei)
         of the jth element bound in molecule i --     -------------- */

  nucl_index = (int **) malloc(atmos.Nmolecule * sizeof(int *));
  for (i = 0;  i < atmos.Nmolecule;  i++) {
    molecule = &atmos.molecules[i];
    nucl_index[i] = (int *) malloc(molecule->Nelement * sizeof(int));
    for (j = 0;  j < molecule->Nelement;  j++) {
      for (nu = 0;  nu < Nnuclei;  nu++) {
	if (nuclei[nu] == &atmos.elements[molecule->pt_index[j]]) {
	  nucl_index[i][j] = nu;
	  break;
	}
      }
    }
  }
  /* --- Number of equations --                        -------------- */

  Nequation = Nnuclei + atmos.Nmolecule;

  /* --- Allocate temporary storage space.
         Quantities for the Newton Raphson --          -------------- */

  f  = (double *) malloc(Nequation * sizeof(double));
  n  = (double *) calloc(Nequation, sizeof(double));
  df = matrix_double(Nequation, Nequation);
  a  = (double *) calloc(Nequation, sizeof(double));

  /* --- Equilibrium constant for each molecule, and neutral
         fraction for each nucleus --                  -------------- */

  Phi = (double *) malloc(atmos.Nmolecule * sizeof(double));
  fn0 = (double *) malloc(Nnuclei * sizeof(double));

  /* --- Initialize structure for Ng convergence acceleration -- ---- */

  Ngn = NgInit(Nequation, Ngdelay=NG_CHEM_DELAY,
	       Ngorder=NG_CHEM_ORDER, Ngperiod=NG_CHEM_PERIOD, n);

  /* --- Go through spatial grid and solve (local) equations -- ----- */

  for (k = 0;  k < atmos.Nspace;  k++) {

    /* --- Collect for each atom the population fraction of
           the neutral stage --                        -------------- */

    for (i = 0;  i < Nnuclei;  i++) {
      if ((atom = nuclei[i]->model) == NULL) {

	/* --- If no atomic model is present for this nucleus -- ---- */

	getfjk(nuclei[i], atmos.ne[k], k, fjk, dfjk);
	fn0[i] = fjk[0];
	a[i]   = nuclei[i]->abund * atmos.nHtot[k];
      } else {

	/* --- Atomic model has been read for this nucleus -- ------- */

	fn0[i] = 0.0;
	for (j = 0;  j < atom->Nlevel;  j++) {
	  if (atom->stage[j] > 0) break;
	  if (atom->active &&
	      atom->initial_solution != OLD_POPULATIONS)
	    fn0[i] += atom->nstar[j][k];
	  else 
	    fn0[i] += atom->n[j][k];
	}
	fn0[i] /= atom->ntotal[k];
	a[i]    = atom->abundance * atmos.nHtot[k];
      }
    }
    PhiHmin = 0.25*pow(CI/atmos.T[k], 1.5) *
      exp(E_ION_HMIN / (KBOLTZMANN * atmos.T[k]));
    fHmin = atmos.ne[k] * fn0[0]*PhiHmin;
    
    /* --- Equilibrium constant for each molecule at this location -- */

    for (i = 0;  i < atmos.Nmolecule;  i++) 
      Phi[i] = equilconstant(&atmos.molecules[i], atmos.T[k]);

    /* --- Initial solution of atomic number densities, and
           the molecules. Assume everything is dissociated -- ------- */

    for (i = 0;  i < Nnuclei;  i++) n[i] = a[i];
    for (i = 0;  i < atmos.Nmolecule;  i++) n[Nnuclei+i] = 0.0;

    /* --- Reset counter in Ng structure, store initial solution - -- */

    Ngn->count = 1;
    for (i = 0;  i < Nequation;  i++)  Ngn->previous[0][i] = n[i];

    /* --- Iterate to convergence, with maximum number NmaxIter -- -- */

    niter = 1;
    while (niter <= NmaxIter) {
      for (i = 0;  i < Nequation;  i++) {
	f[i] = n[i] - a[i];
	for (j = 0;  j < Nequation;  j++) df[i][j] = 0.0;
        df[i][i] = 1.0;
      }
      /* --- Add nHminus to the H number conservation equation -- --- */ 

      f[0] += fHmin * n[0];
      df[0][0] += fHmin;

      /* --- Fill in the rest of the population matrix f[] and its
             derivative df[][] --                      -------------- */

      for (i = 0;  i < atmos.Nmolecule;  i++) {
        molecule = &atmos.molecules[i];
        saha = Phi[i];
	for (j = 0;  j < molecule->Nelement;  j++) {
	  nu = nucl_index[i][j];
          saha *= pow(fn0[nu] * n[nu], molecule->pt_count[j]);

	  /* --- Contributions to equation of conservation for the
                 nuclei in this molecule --            -------------- */

	  f[nu] += molecule->pt_count[j] * n[Nnuclei + i];
	}
	/* --- Saha equation for this molecule --      -------------- */

        saha /= pow(atmos.ne[k], molecule->charge);
        f[Nnuclei + i] -= saha;

        /* --- Fill the derivatives matrix --          -------------- */

        for (j = 0;  j < molecule->Nelement;  j++) {
	  nu = nucl_index[i][j];
          df[nu][Nnuclei + i] += molecule->pt_count[j];
	  df[Nnuclei + i][nu] = -saha * (molecule->pt_count[j]/n[nu]);
	}
      }
      /* --- Solve linearized equations --             -------------- */

      SolveLinearEq(Nequation, df, f, TRUE);
      for (i = 0;  i < Nequation;  i++)  n[i] -= f[i];

      /* --- Check convergence and accelerate if appropriate -- ----- */

      Accelerate(Ngn, n);
      sprintf(messageStr,
	      "\n%s-- Chemical equilibrium: depth %3.3d, iteration %d",
	      (niter == 1) ? "\n" : "", k, niter);

      if ((dnmax = MaxChange(Ngn, messageStr, quiet=TRUE)) <= iterLimit)
	break;
      niter++;
    }
    if (dnmax > iterLimit) {
      sprintf(messageStr, "Iteration not converged:\n"
              " temperature: %6.1f [K], \n"
	      " density: %9.3E [m^-3],\n dnmax: %9.3E\n",
	      atmos.T[k], atmos.nHtot[k], dnmax);
      Error(WARNING, "ChemicalEquilibrium", messageStr);
    }
    /* --- Store population numbers nuclei --          -------------- */

    for (i = 0;  i < Nnuclei;  i++) {
      if ((atom = nuclei[i]->model) != NULL) {
	fraction = n[i] / atom->ntotal[k];

	for (j = 0;  j < atom->Nlevel;  j++) {
	  atom->nstar[j][k] *= fraction;
	  if (atom->n != atom->nstar) atom->n[j][k] *= fraction;
	}
	atom->ntotal[k] = n[i];
      }
    }
    /* --- Store Hmin density --                       -------------- */

    atmos.nHmin[k] = atmos.ne[k] * (n[0] * PhiHmin);

    /* --- Store molecular densities --                -------------- */

    for (i = 0;  i < atmos.Nmolecule;  i++)
      atmos.molecules[i].n[k] = n[Nnuclei + i];
  }

  /* --- Check whether active atom, if present, is in list of nuclei.
         If so print out a warning --                  -------------- */

  for (nu = 0;  nu < Nnuclei;  nu++) {
    atom = nuclei[nu]->model;
    if (atom && atom->active) {
      sprintf(messageStr, "\nReduced number density of"
	      " active atom %s due to molecule%s\n  ",
	      atom->ID, (nuclei[nu]->Nmolecule > 1) ? "s" : "");
      for (i = 0;  i < nuclei[nu]->Nmolecule;  i++) {
	sprintf(tmpStr, "%s%s",
		atmos.molecules[nuclei[nu]->mol_index[i]].ID,
		(i == nuclei[nu]->Nmolecule - 1) ? "\n\n" : ", ");
	strcat(messageStr, tmpStr);
      }
      Error(MESSAGE, routineName, messageStr);
    }
  }

  /* --- Clean up --                                   -------------- */

  free(f);     free(n);     free(a);
  free(Phi);   free(fn0);
  freeMatrix((void **) df);

  NgFree(Ngn);
  free(nuclei);
  for (i = 0;  i < atmos.Nmolecule;  i++) free(nucl_index[i]);
  free(nucl_index);

  if (Nmaxstage) {
    free(fjk);
    free(dfjk);
  }
  getCPU(3, TIME_POLL, "Chemical equilibrium");
}
/* ------- end ---------------------------- ChemicalEquilibrium.c --- */

/* ------- begin -------------------------- partfunction.c ---------- */

double partfunction(struct Molecule *molecule, double T)
{
  register int i;

  double pf = 0.0, t;

  if (T < molecule->Tmin  ||  T > molecule->Tmax)
    return pf;

  /* --- Evaluate polynomial in temperature for partition function -- */

  switch (molecule->fit) {
  case KURUCZ_70:
    pf = molecule->pf_coef[0];
    for (i = 1;  i < molecule->Npf;  i++)
      pf = pf*T + molecule->pf_coef[i];

    pf = exp(pf);
    break;

  case KURUCZ_85:
    t  = T * 1.0E-4;
    pf = molecule->pf_coef[0];
    for (i = 1;  i < molecule->Npf;  i++)
      pf = pf*t + molecule->pf_coef[i];

    pf = exp(pf);
    break;

  case SAUVAL_TATUM_84:
    t = log10(THETA0 / T);
    pf = molecule->pf_coef[0];
    for (i = 1;  i < molecule->Npf;  i++)
      pf = pf*t + molecule->pf_coef[i];

    pf = POW10(pf);
    break;

  case IRWIN_81:
    t = log(T);
    pf = molecule->pf_coef[0];
    for (i = 1;  i < molecule->Npf;  i++)
      pf = pf*t + molecule->pf_coef[i];

    pf = exp(pf);
    break;

  case TSUJI_73:
    break;
  default:
    sprintf(messageStr,
	    "Unknown method for calculation of partition function\n"
            "for molecule %s: %d", molecule->ID, molecule->fit);
    Error(ERROR_LEVEL_2, "partfunction", messageStr);
  }
  return pf;
}
/* ------- end ---------------------------- partfunction.c ---------- */

/* ------- begin -------------------------- equilconstant.c --------- */

double equilconstant(struct Molecule *molecule, double T)
{
  register int i;

  int    mk;
  double kT, t, eqc = 0.0, theta, cgs_to_SI = 1.0;

  if (T < molecule->Tmin  ||  T > molecule->Tmax)
    return eqc;

  /* --- Evaluate polynomial in temperature for equilibrium constant.
         Constant should have units of [m^3,6,9, etc] -- ------------ */

  switch (molecule->fit) {
  case KURUCZ_70:
    kT = KBOLTZMANN * T;
    mk = molecule->Nnuclei - 1 - molecule->charge;
 
    eqc = molecule->eqc_coef[0];
    for (i = 1;  i < molecule->Neqc;  i++)
      eqc = eqc*T + molecule->eqc_coef[i];
    eqc = exp(molecule->Ediss/kT + eqc - 1.5*mk*log(T));

    cgs_to_SI = pow(CUBE(CM_TO_M), mk);
    break;

  case KURUCZ_85:
    t  = T * 1.0E-4;
    kT = KBOLTZMANN * T;
    mk = molecule->Nnuclei - 1 - molecule->charge;

    eqc = molecule->eqc_coef[0];
    for (i = 1;  i < molecule->Neqc;  i++)
      eqc = eqc*t + molecule->eqc_coef[i];
    eqc = exp(molecule->Ediss/kT + eqc - 1.5*mk*log(T));

    cgs_to_SI = pow(CUBE(CM_TO_M), mk);
    break;

  case IRWIN_81:

    /* ---- Use SAUVAL_TATUM_84 for equilibrium constant -- ---------- */
    ;
  case SAUVAL_TATUM_84:
    theta = THETA0 / T;
    t     = log10(theta);
    kT    = KBOLTZMANN * T;

    eqc = molecule->eqc_coef[0];
    for (i = 1;  i < molecule->Neqc;  i++)
      eqc = eqc*t + molecule->eqc_coef[i];
    eqc = POW10((molecule->Ediss/EV) * theta - eqc) * kT;

    /* --- Constants are already given in SI units by Sauval & Tatum - */

    cgs_to_SI = 1.0;
    break;

  case TSUJI_73:
    theta = THETA0 / T;
    kT    = KBOLTZMANN * T;

    eqc = molecule->eqc_coef[0];
    for (i = 1;  i < molecule->Neqc;  i++)
      eqc = eqc*theta + molecule->eqc_coef[i];
    eqc = SQ(kT) * POW10(-eqc);

    cgs_to_SI = pow(CUBE(CM_TO_M) / ERG_TO_JOULE, molecule->Nnuclei-1);
    break;

  default:
    sprintf(messageStr,
	    "Unknown method for calculation of equilibrium constant\n"
            "for molecule %s: %d", molecule->ID, molecule->fit);
    Error(ERROR_LEVEL_2, "equilconstant", messageStr);
  }

  return eqc * cgs_to_SI;
}
/* ------- end ---------------------------- equilconstant.c --------- */
