/* ------- file: -------------------------- formal.c ----------------

       Version:       rh2.0, 1-D spherically symmetric
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Dec 14 16:40:20 2011 --

       --------------------------                      ----------RH-- */

/* --- Formal solution with given source function and PRD emission
       profile. --                                     -------------- */
 
 
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "geometry.h"
#include "error.h"
#include "inputs.h"
#include "background.h"
#include "constant.h"
#include "statistics.h"
#include "xdr.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- Formal.c ---------------- */

double Formal(int nspect, bool_t eval_operator, bool_t redistribute)
{
  const char routineName[] = "Formal";
  register int k, l, mu;

  bool_t   initialize, boundbound, PRD_angle_dep, to_obs, polarized,
           solve_Stokes;
  enum     FeautrierOrder F_order;     
  int      Nspace = atmos.Nspace, Nrays = atmos.Nrays;
  double  *I, *chi, *S, *Psi, *Jdag, dJ, dJmax, *J;
  ActiveSet *as;
  Ray *ray;

  /* --- Retrieve active set as of transitions at wavelength nspect - */

  as = &spectrum.as[nspect];
  alloc_as(nspect, eval_operator);
  
  /* --- Check whether current active set includes a bound-bound
         and/or angledependent PRD transition. Otherwise, only
         angle-independent opacity and source functions are needed -- */ 

  boundbound    = containsBoundBound(as);
  PRD_angle_dep = (containsPRDline(as) && input.PRD_angle_dep);
  polarized     = containsPolarized(as);
  solve_Stokes  = (polarized && input.StokesMode == FULL_STOKES);

  /* --- Allocate temporary space --                   -------------- */

  chi  = (double *) malloc(Nspace * sizeof(double));
  Jdag = (double *) malloc(Nspace * sizeof(double));

  if (eval_operator)
    Psi = (double *) malloc(Nspace * sizeof(double));
  else
    Psi = NULL;

  S = (double *) malloc(Nspace * sizeof(double));
  I = (double *) malloc(Nspace * sizeof(double));

  /* --- Store current mean intensities, initialize new ones to zero- */

  Jdag = (double *) malloc(Nspace * sizeof(double));
  if (input.limit_memory) {
    J = (double *) malloc(atmos.Nspace * sizeof(double));
    readJlambda(nspect, Jdag);
  } else {
    J = spectrum.J[nspect];
    for (k = 0;  k < Nspace;  k++) Jdag[k] = J[k];
  }
  if (spectrum.updateJ)
    for (k = 0;  k < Nspace;  k++) J[k] = 0.0;

  /* --- Case of angle-dependent opacity and source function -- ----- */

  if (polarized || PRD_angle_dep ||
      (atmos.moving && (boundbound || atmos.backgrflags[nspect].hasline))) {
    Error(ERROR_LEVEL_2, routineName,
	  "Angle dependent source functions and opacities\n"
	  " not yet implemented in spherical symmetry");
  } else {

    /* --- The angle-independent case --               -------------- */

    readBackground(nspect, 0, 0);
    Opacity(nspect, 0, 0, initialize=TRUE);
    if (eval_operator) addtoCoupling(nspect);

    for (k = 0;  k < Nspace;  k++) {
      chi[k] = as->chi[k] + as->chi_c[k];
      S[k]   = (as->eta[k] +
		as->eta_c[k] + as->sca_c[k]*Jdag[k]) / chi[k];
    }
    for (mu = 0;  mu < Nrays;  mu++) {
      spectrum.I[nspect][mu] =
	Feautrier(nspect, mu, chi, S, F_order=STANDARD, I, Psi);

      ray = geometry.rays + mu;
      if (eval_operator) {
	for (k = 0;  k < ray->Ns;  k++) Psi[k] /= chi[k];
	addtoGamma_sphere(nspect, ray, I, Psi);
      }
      if (spectrum.updateJ) {
	for (k = 0;  k < ray->Ns;  k++)
	  J[k] += I[k] * ray->wmu[k];
	addtoRates_sphere(nspect, ray, I, redistribute);
      }
    }
  }
  /* --- Write new J for current position in the spectrum -- -------- */

  dJmax = 0.0;
  if (spectrum.updateJ) {
    for (k = 0;  k < Nspace;  k++) {
      dJ = fabs(1.0 - Jdag[k]/J[k]);
      dJmax = MAX(dJmax, dJ);
    }
    if (input.limit_memory) writeJlambda(nspect, J);
  }
  /* --- Clean up --                                 ---------------- */

  if (input.limit_memory) free(J);
  free_as(nspect, eval_operator);

  free(chi);  free(Jdag);
  free(I);    free(S);
  if (eval_operator) free(Psi);

  return dJmax;
}   
/* ------- end ---------------------------- Formal.c ---------------- */

/* ------- begin -------------------------- writeFlux.c ------------- */

bool_t writeFlux(char *fileName)
{
  const char routineName[] = "writeFlux";
  register int mu, nspect;

  bool_t  result = TRUE;
  int     Nrays = geometry.Nrays;
  double *flux, *wmuz, scale_to_surface;
  FILE *fp_flux;
  XDR  xdrs;

  /* --- Write the radiative flux in the z-direction -- ------------ */

  if ((fp_flux = fopen(fileName, "w")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s", fileName);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return FALSE;
  }
  xdrstdio_create(&xdrs, fp_flux, XDR_ENCODE);

  flux = (double *) calloc(spectrum.Nspect, sizeof(double));
  wmuz = (double *) malloc(Nrays * sizeof(double));

  wmuz[0] = 0.25*(SQ(geometry.rays[1].xmu[0]) - SQ(geometry.rays[0].xmu[0]));
  for (mu = 1;  mu < Nrays-1;  mu++) {
    wmuz[mu] = 0.25*(SQ(geometry.rays[mu+1].xmu[0]) -
		     SQ(geometry.rays[mu-1].xmu[0]));
  }
  wmuz[Nrays-1] = 0.25*(SQ(geometry.rays[Nrays-1].xmu[0]) -
			SQ(geometry.rays[Nrays-2].xmu[0]));

  scale_to_surface = SQ((geometry.r[0] + geometry.Radius) / geometry.Radius);
  for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
    for (mu = 0;  mu < geometry.Nrays;  mu++)
      flux[nspect] += spectrum.I[nspect][mu] * wmuz[mu];
    flux[nspect] *= 2.0 * PI * scale_to_surface;
  }
  result &= xdr_vector(&xdrs, (char *) flux, spectrum.Nspect,
		       sizeof(double), (xdrproc_t) xdr_double);
      
  xdr_destroy(&xdrs);
  fclose(fp_flux);
  free(wmuz);  free(flux);

  return result;
}
/* ------- end ---------------------------- writeFlux.c ------------- */

/* ------- begin -------------------------- addtoRates_sphere.c ----- */

void addtoRates_sphere(int nspect, Ray *ray, double *P,
		       bool_t redistribute)
{
  register int n, k, nact;

  int nt;
  double twohnu3_c2, twohc, wmula, *Rij, *Rji;
  ActiveSet *as;
  Atom *atom;
  AtomicLine *line;
  AtomicContinuum *continuum;
  pthread_mutex_t *rate_lock;

  twohc = 2.0*HPLANCK*CLIGHT / CUBE(NM_TO_M);

  as = &spectrum.as[nspect];
  nt = nspect % input.Nthreads;

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
      /* --- Convention: Rij is the rate for transition i -> j -- --- */
 
      if (Rij != NULL) {
	if (input.Nthreads > 1) pthread_mutex_lock(rate_lock);

	for (k = 0;  k < ray->Ns;  k++) {
	  wmula = atom->rhth[nt].Vij[n][k] *
	    atom->rhth[nt].wla[n][k] * ray->wmu[k];
	  Rij[k] += P[k] * wmula;
	  Rji[k] += atom->rhth[nt].gij[n][k] * (twohnu3_c2 + P[k]) * wmula;
	}
	if (input.Nthreads > 1) pthread_mutex_unlock(rate_lock);
      }
    }
  }
}
/* -------- end --------------------------- addtoRates_sphere.c ----- */

/* ------- begin -------------------------- addtoGamma_sphere.c ----- */

void addtoGamma_sphere(int nspect, Ray *ray, double *P, double *Psi)
{
  const char routineName[] = "addtoGamma_sphere";
  register int n, k, m, nact;

  int    i, j, ij, ji, jp, nt;
  double twohnu3_c2, twohc, wmula, *Ieff;

  Atom *atom;
  Molecule *molecule;
  AtomicLine *line;
  AtomicContinuum *continuum;
  MolecularLine *mrt;
  ActiveSet *as;

  twohc = 2.0*HPLANCK*CLIGHT / CUBE(NM_TO_M);

  as = &spectrum.as[nspect];
  nt = nspect % input.Nthreads;

  if (containsActive(as))
    Ieff = (double *) malloc(atmos.Nspace * sizeof(double));

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];

    if (as->Nactiveatomrt[nact] > 0) {
      for (k = 0;  k < atmos.Nspace;  k++)
	Ieff[k] = P[k] - Psi[k] * atom->rhth[nt].eta[k];
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
      /* --- In case of atomic transitions --            -------------- */

      if (input.Nthreads > 1) pthread_mutex_lock(&atom->Gamma_lock);

      ij = i*atom->Nlevel + j;
      ji = j*atom->Nlevel + i;
 
      for (k = 0;  k < ray->Ns;  k++) {
	wmula = atom->rhth[nt].Vij[n][k] *
	  atom->rhth[nt].wla[n][k] * ray->wmu[k];

	atom->Gamma[ji][k] += Ieff[k] * wmula;
	atom->Gamma[ij][k] +=
	  (twohnu3_c2 + Ieff[k]) * atom->rhth[nt].gij[n][k] * wmula;
      }
      /* --- Cross-coupling terms --                   -------------- */
 
      for (k = 0;  k < ray->Ns;  k++) {
        atom->Gamma[ij][k] -= atom->rhth[nt].chi_up[i][k] *
          Psi[k]*atom->rhth[nt].Uji_down[j][k] * ray->wmu[k];
      }
      /* --- If rt->i is also an upper level of another transition that
             is active at this wavelength then Gamma[ij] needs to be
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
          for (k = 0;  k < ray->Ns;  k++) {
            atom->Gamma[ji][k] += atom->rhth[nt].chi_down[j][k] *
              Psi[k]*atom->rhth[nt].Uji_down[i][k] * ray->wmu[k];
          }
        }
      }
      if (input.Nthreads > 1) pthread_mutex_unlock(&atom->Gamma_lock);
    }
  }
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

      for (k = 0;  k < ray->Ns;  k++) {
        if (molecule->n[k]) {
	  wmula = molecule->rhth[nt].Vij[n][k] *
	    molecule->rhth[nt].wla[n][k] * ray->wmu[k];
	  molecule->Gamma[ji][k] += P[k] * wmula;
	  molecule->Gamma[ij][k] += molecule->rhth[nt].gij[n][k] *
	    (twohnu3_c2 + P[k]) * wmula;
	}
      }
      if (input.Nthreads > 1) pthread_mutex_unlock(&molecule->Gamma_lock);
    }
  }

  if (containsActive(as)) free(Ieff);
}
/* ------- end ---------------------------- addtoGamma_sphere.c ----- */

/* ------- begin -------------------------- Hydrostatic.c ----------- */

void Hydrostatic(int NmaxIter, double iterLimit)
{
  const char routineName[] = "Hydrostatic";

  if (atmos.hydrostatic) {
    sprintf(messageStr,
	    "Can only establish hydrostatic equilibrium in 1-D geometry");
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- Hydrostatic.c ----------- */
