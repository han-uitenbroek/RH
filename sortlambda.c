/* ------- file: -------------------------- sortlambda.c ------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Apr 22 09:02:18 2009 --

       --------------------------                      ----------RH-- */

/* --- Sorts wavelengths and determines what transitions are active at
       which wavelength. --                            -------------- */
 
#include <math.h>
#include <stdlib.h>
#include <string.h>


#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "constant.h"
#include "inputs.h"
#include "error.h"
#include "statistics.h"
#include "xdr.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern CommandLine commandline;
extern char messageStr[];


/* ------- begin -------------------------- SortLambda.c ------------ */

void SortLambda()
{
  const char routineName[] = "SortLambda";
  register int kr, n, m, nspect, la, nact;

  bool_t  hunt, unique, result;
  int     Nred, Nspectrum, Nlambda_original, Z, i, j, Nwave;
  double *alpha_original, gbf_0, n_eff, *wavetable;
  ActiveSet *as;
  Atom *atom;
  Molecule *molecule;
  AtomicLine *line;
  AtomicContinuum *continuum;
  MolecularLine *mrt;
  FILE  *fp_wavetable;
  XDR    xdrs;

  getCPU(2, TIME_START, NULL);

  /* --- First read the wavelength table if specified -- ------------ */

  result = TRUE;

  if (strcmp(input.wavetable_input, "none")) {
    if ((fp_wavetable = fopen(input.wavetable_input, "r")) == NULL) {
      sprintf(messageStr, "Unable to open input file %s",
	      input.wavetable_input);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    xdrstdio_create(&xdrs, fp_wavetable, XDR_DECODE);

    result &= xdr_int(&xdrs, &Nwave);
    wavetable = (double *) malloc(Nwave * sizeof(double));
    result &= xdr_vector(&xdrs, (char *) wavetable, Nwave,
			 sizeof(double), (xdrproc_t) xdr_double);
    if (!result) {
      sprintf(messageStr, "Unable to read from input file %s",
	      input.wavetable_input);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    xdr_destroy(&xdrs);
    fclose(fp_wavetable);
  } else
    Nwave = 0;

  /* --- Add reference wavelength if necessary --      -------------- */

  Nspectrum = (atmos.lambda_ref > 0.0) ? Nwave + 1 : Nwave;

  /* --- Go through the active atoms to collect number of wavelengths
         for the lines and continua treated in detail -- ------------ */

  atmos.Nactiveatom = 0;
  for (n = 0;  n < atmos.Natom;  n++) {
    atom = &atmos.atoms[n];
    if (atom->active) {
      for (kr = 0;  kr < atom->Ncont;  kr++)
	Nspectrum += atom->continuum[kr].Nlambda;
      for (kr = 0;  kr < atom->Nline;  kr++)
	Nspectrum += atom->line[kr].Nlambda;

      atom->activeindex = atmos.Nactiveatom;
      atmos.Nactiveatom++;
    }
  }
  /* --- Store the pointers to the active atoms, so that they can be
         enumerated --                                 -------------- */

  if (atmos.Nactiveatom > 0) {
    atmos.activeatoms = (Atom **) malloc(atmos.Nactiveatom *
					 sizeof(Atom *));
    for (n = 0;  n < atmos.Natom;  n++) {
      atom = &atmos.atoms[n];
      if (atom->active)
	atmos.activeatoms[atom->activeindex] = atom;
    }
  } else
    atmos.activeatoms = NULL;

  /* --- Add number of wavelengths for lines in active molecules -- - */

  atmos.Nactivemol = 0;
  for (n = 0;  n < atmos.Nmolecule;  n++) {
    molecule = &atmos.molecules[n];
    if (molecule->active) {
      for (kr = 0;  kr < molecule->Nrt;  kr++)
      Nspectrum += molecule->mrt[kr].Nlambda;

      molecule->activeindex = atmos.Nactivemol;
      atmos.Nactivemol++;
    }
  }
 /* --- Store the pointers to the active molecules, so that they can be
        enumerated --                                  -------------- */

  if (atmos.Nactivemol > 0) {
    atmos.activemols = (Molecule **) malloc(atmos.Nactivemol *
					    sizeof(Molecule *));
    for (n = 0;  n < atmos.Nmolecule;  n++) {
      molecule = &atmos.molecules[n];
      if (molecule->active)
	atmos.activemols[molecule->activeindex] = molecule;
    }
  } else
    atmos.activemols = NULL;

  /* --- Fill the wavelength array --                  -------------- */

  nspect = 0;
  spectrum.lambda = (double *) malloc(Nspectrum * sizeof(double));

  /* --- First the referenece wavelength if specified -- ------------ */

  if (atmos.lambda_ref > 0.0)
    spectrum.lambda[nspect++] = atmos.lambda_ref;

  /* --- Then the wavelength table --                  -------------- */

  for (kr = 0;  kr < Nwave;  kr++)
    spectrum.lambda[nspect++] = wavetable[kr];

  /* --- Finally, all the detailed radiative transitions -- --------- */

  atmos.NPRDactive = 0;

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    for (kr = 0;  kr < atom->Ncont;  kr++) {
      continuum = &atom->continuum[kr];
      for (la = 0;  la < continuum->Nlambda;  la++)
	spectrum.lambda[nspect++] = continuum->lambda[la];
    }
    for (kr = 0;  kr < atom->Nline;  kr++) {
      line = &atom->line[kr];
      for (la = 0;  la < line->Nlambda;  la++)
	spectrum.lambda[nspect++] = line->lambda[la];
      
      if (line->PRD) atmos.NPRDactive++;
    }
  }
  /* --- Active molecular lines --                     -------------- */

  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];
    for (kr = 0;  kr < molecule->Nrt;  kr++) {
      mrt = &molecule->mrt[kr];
      for (la = 0;  la < mrt->Nlambda;  la++)
	spectrum.lambda[nspect++] = mrt->lambda[la];
    }
  }
  /* --- Sort the wavelengths in ascending order --    -------------- */

  qsort(spectrum.lambda, Nspectrum, sizeof(double), qsascend);

  /* --- Check for duplicate wavelengths --            -------------- */

  spectrum.Nspect = 1;
  for (nspect = 1;  nspect < Nspectrum;  nspect++) {
    if (spectrum.lambda[nspect] > spectrum.lambda[nspect-1]) {
      spectrum.lambda[spectrum.Nspect] = spectrum.lambda[nspect];
      spectrum.Nspect++;
    }
  }
  sprintf(messageStr, "\n %s: Found %d unique wavelengths\n",
	  routineName, spectrum.Nspect);
  Error(MESSAGE, routineName, messageStr);
  if (spectrum.Nspect < Nspectrum) {
    sprintf(messageStr, " %s: Eliminated %d duplicate wavelengths\n\n",
	    routineName, Nspectrum - spectrum.Nspect);
    Error(MESSAGE, routineName, messageStr);
  }
  /* --- Allocate space for wavelength array and active sets -- ----- */

  spectrum.lambda = (double *) realloc(spectrum.lambda,
				       spectrum.Nspect*sizeof(double));
  spectrum.as = (ActiveSet *) malloc(spectrum.Nspect * sizeof(ActiveSet));

  /* --- Go through each established wavelength and gather active
         transitions --                                -------------- */

  for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
    as = &spectrum.as[nspect];

    /* --- as->art and as->mrt store the arrays of active
           transitions for each active atom
           and molecule at this wavelength seperately -- ------------ */
 
    if (atmos.Nactiveatom > 0) {
      as->Nactiveatomrt =
	(int *) malloc(atmos.Nactiveatom * sizeof(int));
      as->art = (AtomicTransition **)
	malloc(atmos.Nactiveatom * sizeof(AtomicTransition *));

      for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
	as->Nactiveatomrt[nact] = 0;
	as->art[nact] = (AtomicTransition *)
	  malloc(N_MAX_OVERLAP * sizeof(AtomicTransition));
      }
    } else {
      as->Nactiveatomrt = NULL;
      as->art = NULL;
    }

    if (atmos.Nactivemol > 0) {
      as->Nactivemolrt = (int *)
	malloc(atmos.Nactivemol * sizeof(int));
      as->mrt = (MolTransition **)
	malloc(atmos.Nactivemol * sizeof(MolTransition *));

      for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
	as->Nactivemolrt[nact] = 0;
	as->mrt[nact] = (MolTransition *)
	  malloc(N_MAX_OVERLAP * sizeof(MolTransition));
      }
    } else {
      as->Nactivemolrt = NULL;
      as->mrt = NULL;
    }
  }
  /* --- Determine what transitions are active at which
         wavelengths and store the pointers to those transitions. - - */

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
 
    /*--- Go through the continua first --             -------------- */

    Nred = 0;
    for (kr = 0;  kr < atom->Ncont;  kr++) {
      continuum = &atom->continuum[kr];

      /* --- Store the original wavelength array size -- ------------ */

      Nlambda_original = continuum->Nlambda;

      /* --- Find the indices of the lowest (Nblue) and highest (Nred)
             wavelength of the current transition --     ------------ */

      Hunt(spectrum.Nspect, spectrum.lambda,
	   continuum->lambda[0], &continuum->Nblue);
      Hunt(spectrum.Nspect, spectrum.lambda,
	   continuum->lambda[continuum->Nlambda-1], &Nred);
      continuum->Nlambda = Nred - continuum->Nblue + 1;

      /* --- Store the pointer to the current transition in the
	     active set (as) at each wavelength covered by the current
             transition. Calculate wavelength integration weights - - */

      for (nspect = continuum->Nblue;  nspect <= Nred;  nspect++) {
	as   = &spectrum.as[nspect];
	nact = atom->activeindex;

	as->art[nact][as->Nactiveatomrt[nact]].type = ATOMIC_CONTINUUM;
	as->art[nact][as->Nactiveatomrt[nact]].ptype.continuum = continuum;
	as->Nactiveatomrt[nact]++;
	
	if (as->Nactiveatomrt[nact] == N_MAX_OVERLAP) {
	  sprintf(messageStr,
		  "\n Too many overlapping transitions (> %d) "
		  "for atom %s and nspect = %d\n",
		  as->Nactiveatomrt[nact], atom->ID, nspect);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}
      }
      /* --- In case of Bound-Free transition compute absorption
             cross-section if wavelength dependence is hydrogenic,
	     interpolate if wavelength dependence is given
             explicitly --                             -------------- */
      
      if (continuum->hydrogenic) {
	free(continuum->lambda);
	continuum->lambda = spectrum.lambda + continuum->Nblue;
	continuum->alpha =
	  (double *) realloc(continuum->alpha,
			     continuum->Nlambda*sizeof(double));
	
	Z = atom->stage[continuum->j];
	n_eff = Z * sqrt(E_RYDBERG /
			 (atom->E[continuum->j] - atom->E[continuum->i]));
	gbf_0 = Gaunt_bf(continuum->lambda0, n_eff, Z);
	
	for (la = 0;  la < continuum->Nlambda;  la++) {
	  continuum->alpha[la] = continuum->alpha0 *
	    Gaunt_bf(continuum->lambda[la], n_eff, Z) / gbf_0 *
	    CUBE(continuum->lambda[la]/continuum->lambda0);
	}
      } else {
	alpha_original = continuum->alpha;
	splineCoef(Nlambda_original, continuum->lambda, alpha_original);
	
	continuum->alpha =
	  (double *) malloc(continuum->Nlambda * sizeof(double));
	splineEval(continuum->Nlambda, spectrum.lambda + continuum->Nblue,
		   continuum->alpha, hunt=TRUE);
	
	free(continuum->lambda);
	continuum->lambda = spectrum.lambda + continuum->Nblue;
	free(alpha_original);
      }
    }
    /* --- Then go through the bound-bound transitions -- ----------- */

    for (kr = 0;  kr < atom->Nline;  kr++) {
      line = &atom->line[kr];

      /* --- Store the original wavelength array size --  ----------- */

      Nlambda_original = line->Nlambda;

      /* --- Find the indices of the lowest (Nblue) and highest (Nred)
             wavelength of the current transition --      ----------- */

      Hunt(spectrum.Nspect, spectrum.lambda, line->lambda[0],
	   &line->Nblue);
      Hunt(spectrum.Nspect, spectrum.lambda,
	   line->lambda[line->Nlambda-1], &Nred);
      line->Nlambda = Nred - line->Nblue + 1;

      /* --- Store the pointer to the current transition in the
             active set (as) at each wavelength covered by the current
             transition. Calculate wavelength integration weights - - */

      for (nspect = line->Nblue;  nspect <= Nred;  nspect++) {
	as = &spectrum.as[nspect];
	nact = atom->activeindex;

	as->art[nact][as->Nactiveatomrt[nact]].type = ATOMIC_LINE;
	as->art[nact][as->Nactiveatomrt[nact]].ptype.line = line;
	as->Nactiveatomrt[nact]++;

	if (as->Nactiveatomrt[nact] == N_MAX_OVERLAP) {
	  sprintf(messageStr,
		  "\n Too many overlapping transitions (> %d) "
		  "for atom %s and nspect = %d\n",
		  as->Nactiveatomrt[nact], atom->ID, nspect);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}
      }
      free(line->lambda);
      line->lambda = spectrum.lambda + line->Nblue;
    }
  }

  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];

    Nred = 0;
    for (kr = 0;  kr < molecule->Nrt;  kr++) {
      mrt = &molecule->mrt[kr];

      /* --- Find the indices of the lowest (Nblue) and highest
             wavelength of the current transition and allocate memory
             for the wavelength integration weights -- -------------- */

      Hunt(spectrum.Nspect, spectrum.lambda,
	   mrt->lambda[0], &mrt->Nblue);
      Hunt(spectrum.Nspect, spectrum.lambda,
	   mrt->lambda[mrt->Nlambda-1], &Nred);
      mrt->Nlambda = Nred - mrt->Nblue + 1;
      
      /* --- Repoint to proper position in wavelength array -- ------ */

      free(mrt->lambda);
      mrt->lambda = spectrum.lambda + mrt->Nblue;

      /* --- Store the transition number of the current transition
             in the active molecular set (as) at each wavelength
             covered by the current transition. Calculate wavelength
             integration weights --                      ------------ */

      for (nspect = mrt->Nblue;
	   nspect < mrt->Nblue+mrt->Nlambda;  nspect++) {
	as = &spectrum.as[nspect];
	nact = molecule->activeindex;
	
	as->mrt[nact][as->Nactivemolrt[nact]].type = mrt->type;
	as->mrt[nact][as->Nactivemolrt[nact]].ptype.vrline = mrt;
	as->Nactivemolrt[nact]++;

	if (as->Nactivemolrt[nact] == N_MAX_OVERLAP) {
	  sprintf(messageStr,
		  "\n Too many overlapping transitions (> %d) "
		  "for molecule %s and nspect = %d\n",
		  as->Nactivemolrt[nact], molecule->ID, nspect);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}
      }
    }
  }

  for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
    as = &spectrum.as[nspect];
 
    /* --- For each wavelength in the spectrum gather the unique set of
           lower and upper levels involved in active atomic transitions at that
           wavelength. These are needed in the calculation of the cross
           coupling coefficients in the approximate lambda iteration -- */

    as->Nlower = (int *) malloc(atmos.Nactiveatom * sizeof(int));
    as->Nupper = (int *) malloc(atmos.Nactiveatom * sizeof(int));

    as->upper_levels =
      (int **) malloc(atmos.Nactiveatom * sizeof(int *));
    as->lower_levels =
      (int **) malloc(atmos.Nactiveatom * sizeof(int *));

    for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
      atom = atmos.activeatoms[nact];

      /* --- First, for each wavelength in the spectrum gather
	 the unique set of lower levels of active transitions
	 within each atomic model --               -------------- */

      as->Nlower[nact] = 0;
      as->Nupper[nact] = 0;
      if (as->Nactiveatomrt[nact] > 0) {
	as->lower_levels[nact] =
	  (int *) malloc(as->Nactiveatomrt[nact] * sizeof(int));
	as->upper_levels[nact] =
	  (int *) malloc(as->Nactiveatomrt[nact] * sizeof(int));
      } else {
	as->lower_levels[nact] = NULL;
	as->upper_levels[nact] = NULL;
      }

      for (n = 0;  n < as->Nactiveatomrt[nact];  n++) {
	unique = TRUE;
	switch (as->art[nact][n].type) {
	case ATOMIC_LINE: 
	  i = as->art[nact][n].ptype.line->i;
	  j = as->art[nact][n].ptype.line->j;
	  break;
	case ATOMIC_CONTINUUM:
	  i = as->art[nact][n].ptype.continuum->i;
	  j = as->art[nact][n].ptype.continuum->j;
	  break;
	default:;
	}
	for (m = 0;  m < as->Nlower[nact];  m++) {
	  if (i == as->lower_levels[nact][m]) {
	    unique = FALSE;
	    break;
	  }
	}
	if (unique) {
	  as->lower_levels[nact][as->Nlower[nact]] = i;
	  as->Nlower[nact]++;
	}
 
	/* --- Then add the upper level if unique --   -------------- */

	unique = TRUE;
	for (m = 0;  m < as->Nupper[nact];  m++) {
	  if (j == as->upper_levels[nact][m]) {
	    unique = FALSE;
	    break;
	  }
	}
	if (unique) {
	  as->upper_levels[nact][as->Nupper[nact]] = j;
	  as->Nupper[nact]++;
	}
      }
    }
    /* --- Reallocate space for the atomic transition arrays -- ------- */

    for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
      if (as->Nactiveatomrt[nact] > 0) {
	as->art[nact] = (AtomicTransition *)
	  realloc(as->art[nact],
		  as->Nactiveatomrt[nact] * sizeof(AtomicTransition));

	as->lower_levels[nact] =
	  (int *) realloc(as->lower_levels[nact],
			  as->Nlower[nact] * sizeof(int));
	as->upper_levels[nact] =
	  (int *) realloc(as->upper_levels[nact],
			  as->Nupper[nact] * sizeof(int));
      }
    }
    /* --- Reallocate space for the molecular transition arrays -- -- */

    for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
      if (as->Nactivemolrt[nact] > 0) {
	as->mrt[nact] = (MolTransition *)
	  realloc(as->mrt[nact],
		  as->Nactivemolrt[nact] * sizeof(MolTransition));
      }
    }
  }

  getCPU(2, TIME_POLL, "SortLambda");
}
/* ------- end ---------------------------- SortLambda.c ------------ */
