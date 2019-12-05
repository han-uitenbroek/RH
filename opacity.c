/* ------- file: -------------------------- opacity.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Jan 16 20:02:51 2012 --

       --------------------------                      ----------RH-- */

/* --- Evaluate the total opacity and emissivity due to the active set
       of transitions at the wavelength pointed to by structure as->
       In addition, fill the matrices ->Vij[Nactive][Nspace],
       ->gij[Nactive][Nspace], and ->wla[Nactive][Nspace].

 Note: In moving atmospheres check whether the opacity along the ray
       (to_obs == TRUE) or in the opposite direction
       (to_obs == FALSE) is required.

 Note: If bool initialize is set the gij and Vij of all members of the
       active set are updated. Otherwise, only the bound-bound
       transitions are updated. This prevents unnecessary duplicate
       calculations for bound-free transitions in moving atmospheres.

 Note: If rt->polarized is set the diagonal elements \phi_I of the
       absorption matrix \Phi (see for instance: J.O. Stenflo, 1994,
       in "Solar Magnetic Fields", p. 115.) are used to calculate the
       opacity and emissivity. This allows also a solution in the
       "Polarization-free" approximation.

  See: J. Trujillo Bueno & E. Landi. Degli'Innocenti 1996,
       Solar Physics 164, pp 135-153.
       --                                              -------------- */

 
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "inputs.h"
#include "constant.h"


/* --- Function prototypes --                          -------------- */

double MolProfile(MolecularLine *mrt, int k, int mu, bool_t to_obs,
                  double lambda,
		  double *phi_Q, double *phi_U, double *phi_V,
		  double *psi_Q, double *psi_U, double *psi_V);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- Opacity.c --------------- */

void Opacity(int nspect, int mu, bool_t to_obs, bool_t initialize)
{
  const char routineName[] = "Opacity";
  register int k, n, nact;

  int     la, i, j, vi, vj, lamu, nt, NrecStokes, Nrecphi;
  double *phi, gijk, twohnu3_c2, twohc, hc_4PI, hc_k, hc, fourPI,
    *n_i, *n_j, Bijxhc_4PI, wlambda, chi_l, *chi_Q, *chi_U, *chi_V,
     eta_l, *eta_Q, *eta_U, *eta_V, *chip_Q, *chip_U, *chip_V,
    *phi_Q, *phi_U, *phi_V, *psi_Q, *psi_U, *psi_V;
  bool_t  solveStokes;

  Atom *atom;
  Molecule *molecule;
  AtomicLine *line;
  AtomicContinuum *continuum;
  MolecularLine *mrt;
  ActiveSet *as;

  /* --- Some useful constants --                          ---------- */
 
  hc     = HPLANCK * CLIGHT;
  fourPI = 4.0 * PI;
  hc_4PI = hc / fourPI;
  twohc  = 2.0*hc / CUBE(NM_TO_M);
  hc_k   = hc / (KBOLTZMANN * NM_TO_M);
 
  as = &spectrum.as[nspect];
  nt = nspect % input.Nthreads;

  /* --- If polarized transition is present and we solve for polarized
         radiation we need to fill all four Stokes components -- ---- */

  if (input.StokesMode == FULL_STOKES && containsPolarized(as)) {
    NrecStokes = 4;

    chi_Q = as->chi + atmos.Nspace;
    chi_U = as->chi + 2*atmos.Nspace;
    chi_V = as->chi + 3*atmos.Nspace;

    if (input.magneto_optical) {
      for (k = 0;  k < 3*atmos.Nspace;  k++) as->chip[k] = 0.0;

      chip_Q = as->chip;
      chip_U = as->chip + atmos.Nspace;
      chip_V = as->chip + 2*atmos.Nspace;
    }
  } else
    NrecStokes = 1;

  for (k = 0;  k < NrecStokes*atmos.Nspace;  k++) {
    as->chi[k] = 0.0;
    as->eta[k] = 0.0;
  }

  /* --- Loop over set of active transitions at current wavelength -- */

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];

    /* --- Zero the emissivities for each active atom,
           and the proper thread --                    -------------- */

    if (as->Nactiveatomrt[nact] > 0) {
      for (k = 0;  k < NrecStokes*atmos.Nspace;  k++)
	atom->rhth[nt].eta[k] = 0.0;
    }

    for (n = 0;  n < as->Nactiveatomrt[nact];  n++) {
      switch (as->art[nact][n].type) {
      case ATOMIC_LINE:
	line = as->art[nact][n].ptype.line;
	i = line->i;
	j = line->j;
	n_i = atom->n[i];
	n_j = atom->n[j];

 	/* --- Relative position in line profile --    -------------- */

	la = nspect - line->Nblue;

        /* --- Need all four (or seven) Stokes components -- -------- */

        solveStokes =
	  (line->polarizable && input.StokesMode == FULL_STOKES);

        /* --- Required size of temporary profile array -- ---------- */

	if (input.limit_memory) {
	  if (solveStokes)
	    Nrecphi = (input.magneto_optical) ? 7 : 4;
	  else
	    Nrecphi = 1;

	  phi = (double *) malloc(Nrecphi *
				  atmos.Nspace * sizeof(double));
	}

	if (atmos.moving || solveStokes) {
	  lamu = 2*(atmos.Nrays*la + mu) + to_obs;

	  if (input.limit_memory) {
	    readProfile(line, lamu, phi);
	    if (solveStokes) {
	      phi_Q = phi + atmos.Nspace;
	      phi_U = phi + 2*atmos.Nspace;
	      phi_V = phi + 3*atmos.Nspace;
	      
	      if (input.magneto_optical) {
		psi_Q = phi + 4*atmos.Nspace;
		psi_U = phi + 5*atmos.Nspace;
		psi_V = phi + 6*atmos.Nspace;
	      }
	    }
	  } else {
	    phi = line->phi[lamu];

	    if (solveStokes) {
	      phi_Q = line->phi_Q[lamu];
	      phi_U = line->phi_U[lamu];
	      phi_V = line->phi_V[lamu];

	      if (input.magneto_optical) {
		psi_Q = line->psi_Q[lamu];
		psi_U = line->psi_U[lamu];
		psi_V = line->psi_V[lamu];
	      }
	    }
	  }
	} else {
	  if (input.limit_memory)
	    readProfile(line, la, phi);
	  else
	    phi = line->phi[la];
	}

	twohnu3_c2 = line->Aji / line->Bji;
	gijk = line->Bji / line->Bij;
	Bijxhc_4PI = hc_4PI * line->Bij * line->isotope_frac;
	for (k = 0;  k < atmos.Nspace;  k++) {
	  atom->rhth[nt].gij[n][k] = gijk;
	  atom->rhth[nt].Vij[n][k] = Bijxhc_4PI * phi[k];
	}

	/* --- PRD correction to emission profile --   -------------- */

	if (line->PRD) {
	  if (input.PRD_angle_dep) {
	    lamu = 2*(atmos.Nrays*la + mu) + to_obs;
	    for (k = 0;  k < atmos.Nspace;  k++)
	      atom->rhth[nt].gij[n][k] *= line->rho_prd[lamu][k];
	  } else {
	    for (k = 0;  k < atmos.Nspace;  k++)
	      atom->rhth[nt].gij[n][k] *= line->rho_prd[la][k];
	  }
	}
	/* --- Store wavelength integration weights -- -------------- */

	if (initialize) {
	  wlambda = getwlambda_line(line, la); 
	  for (k = 0;  k < atmos.Nspace;  k++)
	  atom->rhth[nt].wla[n][k] = wlambda * line->wphi[k] / hc_4PI;
	}
	break;

      case ATOMIC_CONTINUUM:
	continuum = as->art[nact][n].ptype.continuum;
	la = nspect - continuum->Nblue;
	i = continuum->i;
	j = continuum->j;
	n_i = atom->n[i];
	n_j = atom->n[j];

	twohnu3_c2 = twohc / CUBE(continuum->lambda[la]);

	/* --- Do not update gij and Vij of bound-free transitions
	       if set has already been initialized --  -------------- */

	if (initialize) {
	  wlambda = getwlambda_cont(continuum, la);
	  for (k = 0;  k < atmos.Nspace;  k++) {
	    atom->rhth[nt].Vij[n][k] = continuum->alpha[la];
	    atom->rhth[nt].gij[n][k] =
	      atom->nstar[i][k] / atom->nstar[j][k] *
	      exp(-hc_k / (continuum->lambda[la] * atmos.T[k]));

	    atom->rhth[nt].wla[n][k] =
	      fourPI/HPLANCK * (wlambda/continuum->lambda[la]);
	  }
	}
	break;

      default:
	sprintf(messageStr, "Invalid transition type");
	Error(ERROR_LEVEL_1, routineName, messageStr);
	twohnu3_c2 = 0.0;
      }
      /* --- Always calculate total opacity and emissivity of set - - */
      
      if (twohnu3_c2) {
	for (k = 0;  k < atmos.Nspace;  k++) {
	  as->chi[k] += atom->rhth[nt].Vij[n][k] *
	    (n_i[k] - atom->rhth[nt].gij[n][k]*n_j[k]);

	  atom->rhth[nt].eta[k] += twohnu3_c2 * atom->rhth[nt].gij[n][k] *
	    atom->rhth[nt].Vij[n][k] * n_j[k];
	}
	/* --- Emission coefficients for Stokes Q, U, V -- ---------- */

	if (as->art[nact][n].type == ATOMIC_LINE  &&  solveStokes) {
	  lamu = 2*(atmos.Nrays*la + mu) + to_obs;

	  eta_Q = atom->rhth[nt].eta + atmos.Nspace;
	  eta_U = atom->rhth[nt].eta + 2*atmos.Nspace;
	  eta_V = atom->rhth[nt].eta + 3*atmos.Nspace;

	  for (k = 0;  k < atmos.Nspace;  k++) {
	    chi_l =
	      Bijxhc_4PI * (n_i[k] - atom->rhth[nt].gij[n][k]*n_j[k]);

	    chi_Q[k] += chi_l * phi_Q[k];
	    chi_U[k] += chi_l * phi_U[k];
	    chi_V[k] += chi_l * phi_V[k];

	    if (input.magneto_optical) {
	      chip_Q[k] += chi_l * psi_Q[k];
	      chip_U[k] += chi_l * psi_U[k];
	      chip_V[k] += chi_l * psi_V[k];
	    }
	    eta_l =
	      Bijxhc_4PI * twohnu3_c2 * atom->rhth[nt].gij[n][k] * n_j[k];

	    eta_Q[k] += eta_l * phi_Q[k];
	    eta_U[k] += eta_l * phi_U[k];
	    eta_V[k] += eta_l * phi_V[k];
	  }
	}
      }
      if (as->art[nact][n].type == ATOMIC_LINE && input.limit_memory)
	free(phi);
    }
  }

  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];

    if (as->Nactivemolrt[nact] > 0) {
      for (k = 0;  k < atmos.Nspace;  k++)
	molecule->rhth[nt].eta[k] = 0.0;
    }

    for (n = 0;  n < as->Nactivemolrt[nact];  n++) {
      switch (as->mrt[nact][n].type) {
      case VIBRATION_ROTATION:
	mrt = as->mrt[nact][n].ptype.vrline;
	vi = mrt->vi;
	vj = mrt->vj;
	n_i = molecule->nv[vi];
	n_j = molecule->nv[vj];

	la = nspect - mrt->Nblue;
	if (atmos.moving) {
	  lamu = 2*(atmos.Nrays*la + mu) + to_obs;
	  phi  = mrt->phi[lamu];
	} else
	  phi = mrt->phi[la];

	if (initialize) {
	  wlambda = getwlambda_mrt(mrt, la);
	  for (k = 0;  k < atmos.Nspace;  k++)
	    molecule->rhth[nt].wla[n][k] = wlambda * mrt->wphi[k] / hc_4PI;
	}

	twohnu3_c2 = mrt->Aji / mrt->Bji;
	Bijxhc_4PI = hc_4PI * mrt->Bij * mrt->gi * mrt->isotope_frac;

	for (k = 0;  k < atmos.Nspace;  k++) {
	  if (molecule->n[k]) {
	    molecule->rhth[nt].Vij[n][k] =
	      Bijxhc_4PI * phi[k] / molecule->pfv[vi][k] *
	      exp(-mrt->Ei / (KBOLTZMANN * atmos.T[k]));
	    molecule->rhth[nt].gij[n][k] =
	      molecule->nvstar[vi][k] / molecule->nvstar[vj][k] *
	      exp(-hc_k / (mrt->lambda0 * atmos.T[k]));
	  }
	}
	break;

      default:
	sprintf(messageStr, "Invalid transition type");
	Error(ERROR_LEVEL_1, routineName, messageStr);
	twohnu3_c2 = 0.0;
      }
      /* --- Always calculate total opacity and emissivity of set - - */
      
      if (twohnu3_c2) {
	for (k = 0;  k < atmos.Nspace;  k++) {
	  as->chi[k] += molecule->rhth[nt].Vij[n][k] *
	    (n_i[k] - molecule->rhth[nt].gij[n][k]*n_j[k]);

	  molecule->rhth[nt].eta[k] +=
	    twohnu3_c2 * molecule->rhth[nt].gij[n][k] *
	    molecule->rhth[nt].Vij[n][k] * n_j[k];
	}
      }
    }
  }

  /* --- Add all the active contributions into the total
         emissivities of this active set --            -------------- */

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    if (as->Nactiveatomrt[nact] > 0) {
      for (k = 0;  k < NrecStokes*atmos.Nspace;  k++)
        as->eta[k] += atom->rhth[nt].eta[k];
    }
  }
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) { 
    molecule = atmos.activemols[nact];
    if (as->Nactivemolrt[nact] > 0) {
      for (k = 0;  k < atmos.Nspace;  k++)
        as->eta[k] += molecule->rhth[nt].eta[k];
    }
  }
}
/* ------- end ---------------------------- Opacity.c --------------- */

/* ------- begin -------------------------- alloc_as.c -------------- */

void alloc_as(int nspect, bool_t crosscoupling)
{
  register int n, m, nact;

  int i, j, NrecStokes, NrecStokes_as, Nactive, nt;
  Atom *atom;
  Molecule *molecule;
  ActiveSet *as;

  as = &spectrum.as[nspect];
  nt = nspect % input.Nthreads;

  /* --- Allocate space for background opacities and emissivity -- -- */

  if (atmos.backgrflags[nspect].ispolarized &&
      input.StokesMode == FULL_STOKES) {
    NrecStokes = 4;

    if (input.magneto_optical)
      as->chip_c = (double *) malloc(3*atmos.Nspace * sizeof(double));
  } else
    NrecStokes = 1;

  as->chi_c = (double *) malloc(NrecStokes*atmos.Nspace * sizeof(double));
  as->eta_c = (double *) malloc(NrecStokes*atmos.Nspace * sizeof(double));
  as->sca_c = (double *) malloc(atmos.Nspace * sizeof(double));

  /* --- Now for the active part --                    -------------- */

  if (input.StokesMode == FULL_STOKES  && containsPolarized(as)) {
    NrecStokes_as = 4;

    if (input.magneto_optical)
      as->chip = (double *) malloc(3*atmos.Nspace * sizeof(double));
  } else
    NrecStokes_as = 1;

  as->chi = (double *) malloc(NrecStokes_as*atmos.Nspace * sizeof(double));
  as->eta = (double *) malloc(NrecStokes_as*atmos.Nspace * sizeof(double));

  /* --- Allocate memory for the emissivity in each active
         atom seperately. The contributions are needed individually
         in the approximate operator --                -------------- */

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];

    if (as->Nactiveatomrt[nact] > 0) {
      atom->rhth[nt].eta =
	(double *) malloc(NrecStokes_as * atmos.Nspace * sizeof(double));

      atom->rhth[nt].Vij =
	matrix_double(as->Nactiveatomrt[nact], atmos.Nspace);
      atom->rhth[nt].gij =
	matrix_double(as->Nactiveatomrt[nact], atmos.Nspace);
      atom->rhth[nt].wla =
	matrix_double(as->Nactiveatomrt[nact], atmos.Nspace);

      /* --- Allocate pointer space for cross-coupling coefficients - */

      if (crosscoupling) {
	atom->rhth[nt].chi_down = 
	  (double **) calloc(atom->Nlevel, sizeof(double *));
	atom->rhth[nt].chi_up   =
	  (double **) calloc(atom->Nlevel, sizeof(double *));
	atom->rhth[nt].Uji_down =
	  (double **) calloc(atom->Nlevel, sizeof(double *));

	for (m = 0;  m < as->Nlower[nact];  m++) {
	  i = as->lower_levels[nact][m];
	  atom->rhth[nt].chi_up[i] =
	    (double *) malloc(atmos.Nspace * sizeof(double));
	}
	for (m = 0;  m < as->Nupper[nact];  m++) {
	  j = as->upper_levels[nact][m];
	  atom->rhth[nt].chi_down[j] =
	    (double *) malloc(atmos.Nspace * sizeof(double));
	  atom->rhth[nt].Uji_down[j] =
	    (double *) malloc(atmos.Nspace * sizeof(double));
	}
      }
    }
  }
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];

    if (as->Nactivemolrt[nact] > 0) {
      molecule->rhth[nt].eta =
	(double *) malloc(atmos.Nspace * sizeof(double));

      molecule->rhth[nt].Vij = 
	matrix_double(as->Nactivemolrt[nact], atmos.Nspace);
      molecule->rhth[nt].gij =
	matrix_double(as->Nactivemolrt[nact], atmos.Nspace);
      molecule->rhth[nt].wla =
	matrix_double(as->Nactivemolrt[nact], atmos.Nspace);
    }
  }
}
/* ------- end ---------------------------- alloc_as.c -------------- */

/* ------- begin -------------------------- free_as.c --------------- */

void free_as(int nspect, bool_t crosscoupling)
{
  register int nact, n, m;

  int i, j, nt;
  Atom *atom;
  Molecule *molecule;
  ActiveSet *as;

  as = &spectrum.as[nspect];
  nt = nspect % input.Nthreads;

  free(as->chi_c);
  free(as->eta_c);
  free(as->sca_c);

  free(as->chi);
  free(as->eta);

  if (input.StokesMode == FULL_STOKES && input.magneto_optical) {
    if (atmos.backgrflags[nspect].ispolarized)
      free(as->chip_c);
    if (containsPolarized(as))
      free(as->chip);
  }

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];

    if (as->Nactiveatomrt[nact] > 0) {
      free(atom->rhth[nt].eta);

      freeMatrix((void **) atom->rhth[nt].Vij);
      freeMatrix((void **) atom->rhth[nt].gij);
      freeMatrix((void **) atom->rhth[nt].wla);

      if (crosscoupling) {
	for (m = 0;  m < as->Nlower[nact];  m++) {
	  i = as->lower_levels[nact][m];
	  free(atom->rhth[nt].chi_up[i]);
	  atom->rhth[nt].chi_up[i] = NULL;
	}
	free(atom->rhth[nt].chi_up);

	for (m = 0;  m < as->Nupper[nact];  m++) {
	  j = as->upper_levels[nact][m];
	  free(atom->rhth[nt].chi_down[j]);
	  free(atom->rhth[nt].Uji_down[j]);

	  atom->rhth[nt].chi_down[j] = NULL;
	  atom->rhth[nt].Uji_down[j] = NULL;
	}
	free(atom->rhth[nt].chi_down);
	free(atom->rhth[nt].Uji_down);
      }
    }
  }
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];
    if (as->Nactivemolrt[nact] > 0) {
      free(molecule->rhth[nt].eta);

      freeMatrix((void **) molecule->rhth[nt].Vij);
      freeMatrix((void **) molecule->rhth[nt].gij);
      freeMatrix((void **) molecule->rhth[nt].wla);
    }
  }
}
/* ------- end ---------------------------- free_as.c --------------- */

/* ------- begin -------------------------- containsPolarized.c ----- */

bool_t containsPolarized(ActiveSet *as)
{
  register int n, nact;

  if (!atmos.Stokes || input.StokesMode == FIELD_FREE) return FALSE;

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    for (n = 0;  n < as->Nactiveatomrt[nact];  n++) {
      if (as->art[nact][n].type == ATOMIC_LINE &&
	  as->art[nact][n].ptype.line->polarizable) {
	return TRUE;
      }
    }
  }
  return FALSE;
}
/* ------- end ---------------------------- containsPolarized.c ----- */

/* ------- begin -------------------------- containsBoundBound.c ---- */

bool_t containsBoundBound(ActiveSet *as)
{
  register int n, nact;

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    for (n = 0;  n < as->Nactiveatomrt[nact];  n++) {
      if (as->art[nact][n].type == ATOMIC_LINE) {
	return TRUE;
      }
    }
  }
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    for (n = 0;  n < as->Nactivemolrt[nact];  n++) {
      if (as->mrt[nact][n].type == VIBRATION_ROTATION) {
	return TRUE;
      }
    }
  }
  return FALSE;
}
/* ------- end ---------------------------- containsBoundBound.c ---- */

/* ------- begin -------------------------- containsActive.c -------- */

bool_t containsActive(ActiveSet *as)
{
  register int n, nact;

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    if (as->Nactiveatomrt[nact] > 0)
	return TRUE;
  }

  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    if (as->Nactivemolrt[nact] > 0)
	return TRUE;
  }

  return FALSE;
}
/* ------- end ---------------------------- containsActive.c -------- */

/* ------- begin -------------------------- containsPRDline.c ------- */

bool_t containsPRDline(ActiveSet *as)
{
  register int n, nact;

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    for (n = 0;  n < as->Nactiveatomrt[nact];  n++) {
      if (as->art[nact][n].type == ATOMIC_LINE &&
	  as->art[nact][n].ptype.line->PRD) {
	return TRUE;
      }
    }
  }
  return FALSE;
}
/* ------- end ---------------------------- containsPRDline.c ------- */

/* ------- begin -------------------------- mrt_locate.c ------------ */
 
void mrt_locate(int N, MolecularLine *lines, double lambda, int *low)
{
  int high, index, increment;

  if ((*low <= 0)  ||  (*low > N-1)) {

    /* --- Input guess not useful here, go to bisection --  --------- */

    *low = 0;
    high = N;
  } else {

    /* --- Else hunt up or down to bracket value --    -------------- */ 

    increment = 1;
    if (lambda >= lines[*low].lambda0) {
      high = *low + increment;
      if (*low == N-1) return;

      /* --- Hunt up --                                -------------- */

      while (lambda >= lines[high].lambda0) {
	*low = high;
	increment += increment;
	high = *low + increment;
        if (high >= N) { high = N;  break; }
      }
    } else {
      high = *low;
      if (*low == 0) return;

      /* --- Hunt down --                              -------------- */

      while (lambda <= lines[*low].lambda0) {
	high = *low;
	increment += increment;
	*low = high - increment;
        if (*low <= 0) { *low = 0;  break; }
      }
    }
  }
  /* --- Bisection algorithm --                        -------------- */

  while (high - *low > 1) {
    index = (high + *low) >> 1;
    if (lambda >= lines[index].lambda0)
      *low = index;
    else
      high = index;
  }
}
/* ------- end ---------------------------- mrt_locate.c ------------ */

/* ------- begin -------------------------- MolecularOpacity.c ------ */

flags MolecularOpacity(double lambda, int nspect, int mu, bool_t to_obs,
		       double *chi, double *eta, double *chip)
{
  register int n, k, kr;

  /* --- Opacity due to molecular lines in the background. LTE
         populations are assumed. If magnetic fields are present
         the Stokes emission coefficients eta_{Q,U,V} are also
         calculated.
         --                                            -------------- */

  int    NrecStokes;
  double dlamb_char0, dlamb_charN, kT, hc_la, hc, fourPI, hc_4PI,
    Bijhc_4PI, ni_gi, nj_gj, twohnu3_c2, phi, phi_Q, phi_U, phi_V,
    psi_Q, psi_U, psi_V, chi_l, *chi_Q, *chi_U, *chi_V,
    eta_l, *eta_Q, *eta_U, *eta_V, *chip_Q, *chip_U, *chip_V;
  Molecule        *molecule;
  MolecularLine   *mrt;
  flags            backgrflags;

  backgrflags.hasline     = FALSE;
  backgrflags.ispolarized = FALSE;

  hc     = HPLANCK * CLIGHT;
  fourPI = 4.0 * PI;
  hc_4PI = hc / fourPI;

  /* --- Initialize the contributions for this wavelength/angle -- -- */

  if (atmos.Stokes) {
    NrecStokes = 4;

    /* --- Use pointers to sub-arrays for Q, U, and V -- ------------ */

    chi_Q = chi + atmos.Nspace;
    chi_U = chi + 2*atmos.Nspace;
    chi_V = chi + 3*atmos.Nspace;

    eta_Q = eta + atmos.Nspace;
    eta_U = eta + 2*atmos.Nspace;
    eta_V = eta + 3*atmos.Nspace;

    if (input.magneto_optical) {
      chip_Q = chip;
      chip_U = chip + atmos.Nspace;
      chip_V = chip + 2*atmos.Nspace;

      for (k = 0;  k < 3*atmos.Nspace;  k++) chip[k] = 0.0;
    }
  } else
    NrecStokes = 1;

  for (k = 0;  k < NrecStokes*atmos.Nspace;  k++) {
    chi[k] = 0.0;
    eta[k] = 0.0;
  }

  /* --- Loop through the list of molecules and their lines -- ------ */

  for (n = 0;  n < atmos.Nmolecule;  n++) {
    molecule = &atmos.molecules[n];

    if ((molecule->Nrt > 0  &&  !molecule->active)) {
      dlamb_char0 = lambda * molecule->mrt[0].qwing *
	(atmos.vmicro_char / CLIGHT);
      dlamb_charN = lambda * molecule->mrt[molecule->Nrt-1].qwing *
	(atmos.vmicro_char / CLIGHT);

      if (lambda >= molecule->mrt[0].lambda0 - dlamb_char0 &&
	  lambda <= molecule->mrt[molecule->Nrt-1].lambda0 + dlamb_charN) {

	for (kr = 0;  kr < molecule->Nrt;  kr++) {
	  mrt = &molecule->mrt[kr];
          dlamb_char0 = lambda * mrt->qwing * (atmos.vmicro_char / CLIGHT);

	  if (fabs(mrt->lambda0 - lambda) <= dlamb_char0) {
	    hc_la      = (HPLANCK * CLIGHT) / (mrt->lambda0 * NM_TO_M);
	    Bijhc_4PI  = hc_4PI * mrt->Bij * mrt->isotope_frac * mrt->gi;
	    twohnu3_c2 = mrt->Aji / mrt->Bji;

	    backgrflags.hasline = TRUE;
	    if (mrt->polarizable) {
	      backgrflags.ispolarized = TRUE;
	      if (mrt->zm == NULL) mrt->zm = MolZeeman(mrt);
	    }

	    for (k = 0;  k < atmos.Nspace;  k++) {
	      if (molecule->n[k] > 0.0) {
                phi = MolProfile(mrt, k, mu, to_obs, lambda,
				 &phi_Q, &phi_U, &phi_V,
				 &psi_Q, &psi_U, &psi_V);

		kT    = 1.0 / (KBOLTZMANN * atmos.T[k]);
		ni_gi = molecule->n[k] * exp(-mrt->Ei * kT) /
		  molecule->pf[k];
                nj_gj = ni_gi * exp(-hc_la * kT);

                chi_l = Bijhc_4PI * (ni_gi - nj_gj);
		eta_l = Bijhc_4PI * twohnu3_c2 * nj_gj;

		chi[k] += chi_l * phi;
		eta[k] += eta_l * phi;

		if (mrt->zm != NULL) {
		  chi_Q[k] += chi_l * phi_Q;
		  chi_U[k] += chi_l * phi_U;
		  chi_V[k] += chi_l * phi_V;

		  eta_Q[k] += eta_l * phi_Q;
		  eta_U[k] += eta_l * phi_U;
		  eta_V[k] += eta_l * phi_V;

		  if (input.magneto_optical) {
		    chip_Q[k] += chi_l * psi_Q;
		    chip_U[k] += chi_l * psi_U;
		    chip_V[k] += chi_l * psi_V;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return backgrflags;
}
/* ------- end ---------------------------- MolecularOpacity.c ------ */

/* ------- begin -------------------------- MolProfile.c ------------ */

double MolProfile(MolecularLine *mrt, int k, int mu, bool_t to_obs,
                  double lambda,
		  double *phi_Q, double *phi_U, double *phi_V,
		  double *psi_Q, double *psi_U, double *psi_V)
{
  register int nz;

  double    v, phi_sm, phi_sp, phi_pi, psi_sm, psi_sp, psi_pi, adamp,
            vB, H, F, sv, phi_sigma, phi_delta, sign, sin2_gamma, phi,
            psi_sigma, psi_delta;
  Molecule *molecule = mrt->molecule;

  /* --- Returns the normalized profile for a molecular line,
         and calculates the Stokes profile components if necessary -- */

  adamp = mrt->Aji * (mrt->lambda0 * NM_TO_M) / (4.0*PI *
						 molecule->vbroad[k]);
  v = (lambda/mrt->lambda0 - 1.0) * CLIGHT/molecule->vbroad[k];
  if (atmos.moving) {
    if (to_obs)
      v += vproject(k, mu) / molecule->vbroad[k];
    else
      v -= vproject(k, mu) / molecule->vbroad[k];
  }
  sv = 1.0 / (SQRTPI * molecule->vbroad[k]);

  if (mrt->polarizable) {
    sin2_gamma = 1.0 - SQ(atmos.cos_gamma[mu][k]);
    vB   = (LARMOR * mrt->lambda0) * atmos.B[k] / molecule->vbroad[k];
    sign = (to_obs) ? 1.0 : -1.0;

    phi_sm = phi_pi = phi_sp = 0.0;
    psi_sm = psi_pi = psi_sp = 0.0;

    for (nz = 0;  nz < mrt->zm->Ncomponent;  nz++) {
      H = Voigt(adamp, v - mrt->zm->shift[nz]*vB, &F, HUMLICEK);

      switch (mrt->zm->q[nz]) {
      case -1:
	phi_sm += mrt->zm->strength[nz] * H;
	psi_sm += mrt->zm->strength[nz] * F;
	break;
      case  0:
	phi_pi += mrt->zm->strength[nz] * H;
	psi_pi += mrt->zm->strength[nz] * F;
	break;
      case  1:
	phi_sp += mrt->zm->strength[nz] * H;
	psi_sp += mrt->zm->strength[nz] * F;
      }
    }
    phi_sigma = phi_sp + phi_sm;
    phi_delta = 0.5*phi_pi - 0.25*phi_sigma;

    phi = (phi_delta*sin2_gamma + 0.5*phi_sigma) * sv;

    *phi_Q = sign * phi_delta * sin2_gamma * atmos.cos_2chi[mu][k] * sv;
    *phi_U = phi_delta * sin2_gamma * atmos.sin_2chi[mu][k] * sv;
    *phi_V = sign * 0.5*(phi_sp - phi_sm) * atmos.cos_gamma[mu][k] * sv;

    if (input.magneto_optical) {
      psi_sigma = psi_sp + psi_sm;
      psi_delta = 0.5*psi_pi - 0.25*psi_sigma;

      *psi_Q = sign * psi_delta * sin2_gamma * atmos.cos_2chi[mu][k] * sv;
      *psi_U = psi_delta * sin2_gamma * atmos.sin_2chi[mu][k] * sv;
      *psi_V = sign * 0.5*(psi_sp - psi_sm) * atmos.cos_gamma[mu][k] * sv;
    }
  } else
   phi = Voigt(adamp, v, NULL, ARMSTRONG) * sv;

  return phi;
}
/* ------- end ---------------------------- MolProfile.c ------------ */
