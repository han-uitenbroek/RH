/* ------- file: -------------------------- readatom.c --------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Feb 26 16:11:25 2014 --

       --------------------------                      ----------RH-- */

/* --- Reads atomic model into Atom structure pointed to by atom.

       Set keyword active == FALSE if the atomic model is to be used for
       background opacity calculations. In this case no attempt will be
       made to read collisional data.

       Even for passive atoms NLTE populations can be read.

 Note: Element abundance is taken from the abundance input file rather
       than the atomic data file.

 Note: - For LTE metals metal.n is just an alias for metal.nstar.
       - For hydrogen memory for H.n gets allocated in distribute_nH
       (hydrogen.c) when atmos.H_LTE is false and atmos.H is not active.

       --                                              -------------- */

 
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "statistics.h"

#define  COMMENT_CHAR   "#"
#define  MAX_ABUND_ERROR    0.001

/* --- Function prototypes --                          -------------- */

void  distribute_nH(void);
char *getAtomID(char *atom_file);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern CommandLine commandline;
extern char messageStr[];


/* ------- begin -------------------------- readAtom.c -------------- */

void readAtom(Atom *atom, char *atom_file, bool_t active)
{
  const char routineName[] = "readAtom";
  register int kr, krp, kf, la, k, n;

  char    inputLine[MAX_LINE_SIZE], shapeStr[20], vdWstr[20], nuDepStr[20],
          symmStr[20], optionStr[20], labelStr[MAX_LINE_SIZE];
  bool_t  Debeye, exit_on_EOF, match;
  int     i, j, Nlevel, Nrad, Nline, Ncont, Nfixed, 
          Nspace = atmos.Nspace,
          Nread, Nrequired, checkPoint, L, nq, status;
  double  f, C, lambda0, lambdamin, vtherm, S, Ju, Jl,
    c_sum, waveratio, lambda_air;

  AtomicLine *line, *line1;
  AtomicContinuum *continuum;
  FixedTransition *fixed;

  getCPU(3, TIME_START, NULL);

  C = 2*PI * (Q_ELECTRON/EPSILON_0) * (Q_ELECTRON/M_ELECTRON) / CLIGHT;

  /* --- Open the data file for current model atom --  -------------- */

  initAtom(atom);
  if ((atom->fp_input = fopen(atom_file, "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", atom_file);
    Error(ERROR_LEVEL_2, routineName, atom_file);
  } else {
    sprintf(messageStr, " -- reading input file: %s %s",
	    atom_file, (active) ? "(active)\n\n" : "(passive)\n");
    Error(MESSAGE, routineName, messageStr);
  }
  atom->active = active;

  /* --- Read atom ID and convert to uppercase --     -------------- */
 
  getLine(atom->fp_input, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(inputLine, "%2s", atom->ID);
  checkNread(Nread, Nrequired=1, routineName, checkPoint=1);
  for (n = 0;  n < (int) strlen(atom->ID);  n++)
    atom->ID[n] = toupper(atom->ID[n]);
  if (strlen(atom->ID) == 1) strcat(atom->ID, " ");
 
  /* --- NOTE: atomic weight and abundance are read from the 
         abundance input file (abundance.input by default).
 
         When atom is part of the background the abundance and atomic
         weight are taken from the list of elements in atmos.
         See: abundance.c --                          -------------- */

  match = FALSE;
  for (n = 0;  n < atmos.Nelem;  n++) {
    if (strstr(atmos.elements[n].ID, atom->ID)) {
      if (atmos.elements[n].abundance_set) {
	atom->periodic_table = n;
	atom->abundance = atmos.elements[n].abund;
	atom->weight    = atmos.elements[n].weight;
	match = TRUE;
      }
      break;
    }
  }
  if (!match) {
    sprintf(messageStr, " No matching element in periodic table for "
	    " element %s in file %s, or abundance not specified",
	    atom->ID, atom_file);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  /* --- Get Number of levels, lines fixed transitions, and continua  */
 
  getLine(atom->fp_input, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(inputLine, "%d %d %d %d",
		 &atom->Nlevel, &atom->Nline, &atom->Ncont, &atom->Nfixed);
  checkNread(Nread, Nrequired=4, routineName, checkPoint=2);
  Nlevel = atom->Nlevel;
  Nline  = atom->Nline;  Ncont  = atom->Ncont;  Nrad = Nline + Ncont;
  Nfixed = atom->Nfixed;

  atom->E = (double *) malloc(Nlevel * sizeof(double));
  atom->g = (double *) malloc(Nlevel * sizeof(double));
  atom->stage = (int *)   malloc(Nlevel * sizeof(int));
  atom->label = (char **) malloc(Nlevel * sizeof(char *));

  /* --- Read in the level energies, statistical weights, labels,
         and ionization stage --                       -------------- */

  for (i = 0;  i < Nlevel;  i++) {
    atom->label[i] = (char *) calloc((ATOM_LABEL_WIDTH+1), sizeof(char)); 
    getLine(atom->fp_input, COMMENT_CHAR, inputLine , exit_on_EOF=TRUE);
    Nread = sscanf(inputLine, "%lf %lf '%20c' %d",
      &atom->E[i], &atom->g[i], atom->label[i], &atom->stage[i]);
    checkNread(Nread, Nrequired=4, routineName, checkPoint=3);

    atom->E[i] *= (HPLANCK * CLIGHT) / CM_TO_M;
  }
  if (atom->stage[Nlevel-1] != (atom->stage[Nlevel-2] + 1)) {
    sprintf(messageStr,
	    "Atomic model %s in file %s does not have overlying continuum",
	    atom->ID, atom_file);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  atom->nstar  = matrix_double(Nlevel, Nspace);
  atom->ntotal = (double *) malloc(Nspace * sizeof(double));

  for (k = 0;  k < Nspace;  k++)
    atom->ntotal[k] = atom->abundance * atmos.nHtot[k];

  /* --- Ratio of thermal velocity and speed of light for use in
         Doppler width for this particular atomic weight --  -------- */

  if (atom->Nline > 0) {
    atom->vbroad = (double *) malloc(Nspace * sizeof(double));
    vtherm = 2.0*KBOLTZMANN/(AMU * atom->weight);
    for (k = 0;  k < Nspace;  k++)
      atom->vbroad[k] = sqrt(vtherm*atmos.T[k] + SQ(atmos.vturb[k]));
  }

  /* --- Check validity of input.isum for active atom -- ------------ */

  if (atom->active  &&  (input.isum < -1  ||  input.isum >= Nlevel)) {
    sprintf(messageStr, "Invalid value for ISUM = %d, Nlevel = %d",
	    input.isum, Nlevel);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  /* --- Go through the bound-bound transitions --     -------------- */

  atom->Nprd = 0;
  atom->line = (AtomicLine *) malloc(Nline * sizeof(AtomicLine));

  for (kr = 0;  kr < Nline;  kr++) {
    line = atom->line + kr;
    initAtomicLine(line);
    line->atom = atom;
    line->isotope_frac = 1.0;

    getLine(atom->fp_input, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
    Nread = sscanf(inputLine,
		   "%d %d %lf %s %d %s %lf %lf %s %lf %lf %lf %lf %lf %lf %lf",
		   &j, &i, &f, shapeStr, &line->Nlambda, symmStr,
                   &line->qcore, &line->qwing, vdWstr,
		   &line->cvdWaals[0], &line->cvdWaals[1],
		   &line->cvdWaals[2], &line->cvdWaals[3],
		   &line->Grad, &line->cStark, &line->g_Lande_eff);
    checkNread(Nread, Nrequired=15, routineName, checkPoint=4);
    if (Nread == 15) line->g_Lande_eff = 0.0;

    line->j = MAX(i, j);  line->i = MIN(i, j);
    i = line->i;
    j = line->j;

    lambda0 = (HPLANCK * CLIGHT) / (atom->E[j] - atom->E[i]);
    line->Aji = C / SQ(lambda0) * (atom->g[i] / atom->g[j]) * f;
    line->Bji = CUBE(lambda0) / (2.0 * HPLANCK * CLIGHT) * line->Aji;
    line->Bij = (atom->g[j] / atom->g[i]) * line->Bji;
    line->lambda0 = lambda0 / NM_TO_M;

    /* --- Parse the shape string --                   -------------- */

    if (!strstr(shapeStr, "PRD") &&
	!strstr(shapeStr, "VOIGT") &&
	!strstr(shapeStr, "GAUSS") &&
	!strstr(shapeStr, "COMPOSIT")) {
      sprintf(messageStr, "Invalid value for line-shape string: %s",
	      shapeStr);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    if (strstr(shapeStr, "PRD") && input.PRD_NmaxIter > 0) {
      atom->Nprd++;
      line->PRD = TRUE;

      if (input.limit_memory && atom->active) {
        sprintf(messageStr,
                "Cannot limit memory with PRD lines present in active atom");
        Error(ERROR_LEVEL_2, routineName, messageStr);
      }
    }

    if (strstr(shapeStr, "GAUSS")) line->Voigt = FALSE;

    if (strstr(shapeStr, "COMPOSIT")) {
      getLine(atom->fp_input, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
      Nread = sscanf(inputLine, "%d", &line->Ncomponent);
      line->c_shift = (double *) malloc(line->Ncomponent * sizeof(double));
      line->c_fraction =
	(double *) malloc(line->Ncomponent * sizeof(double));

      c_sum = 0.0;
      for (n = 0;  n < line->Ncomponent;  n++) {
	getLine(atom->fp_input, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
	Nread = sscanf(inputLine, "%lf %lf",
		       &line->c_shift[n], &line->c_fraction[n]);
	c_sum += line->c_fraction[n];
      }
      if (c_sum > (1.0 + MAX_ABUND_ERROR) ||
	  c_sum < (1.0 - MAX_ABUND_ERROR)) {
        vacuum_to_air(1, &line->lambda0, &lambda_air);
	sprintf(messageStr,
		"Line %d -> %d (%9.3f [nm]): \n"
		" Component fractions do not add up to unity: %f",
		line->j, line->i, lambda_air, c_sum);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      } else {
	vacuum_to_air(1, &line->lambda0, &lambda_air);
	sprintf(messageStr,
		" Line %d -> %d (%9.3f [nm]) has %d components\n",
		line->j, line->i, lambda_air, line->Ncomponent);
	Error(MESSAGE, routineName, messageStr);
      }
    } else {
      line->Ncomponent = 1;

      line->c_shift = (double *) malloc(sizeof(double));
      line->c_shift[0] = 0.0;
      line->c_fraction = (double *) malloc(sizeof(double));
      line->c_fraction[0] = 1.0;
    }

    if (strstr(vdWstr, "PARAMTR"))
      line->vdWaals = RIDDER_RENSBERGEN;
    else if (strstr(vdWstr, "UNSOLD")) {
      line->vdWaals = UNSOLD;
      line->cvdWaals[3] = line->cvdWaals[1] = 0.0;
    } else if (strstr(vdWstr, "BARKLEM")) {
      line->vdWaals = BARKLEM;
      if (!getBarklemactivecross(line)) {
	sprintf(messageStr,
		"Line %3d -> %3d: cannot treat line "
		"with Barklem type broadening. Using UNSOLD.", j, i);
	Error(WARNING, routineName, messageStr);
	line->vdWaals = UNSOLD;
	line->cvdWaals[3] = line->cvdWaals[1] = 0.0;
      }
    } else {
      sprintf(messageStr, "Invalid value for vd Waals string: %s", vdWstr);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
 
    line->symmetric   = (strstr(symmStr, "ASYMM")) ? FALSE : TRUE;
    line->polarizable = FALSE;

    if (atom->active) {

      /* --- Allocate space for up- and downward radiative rates -- - */

      line->Rij = (double *) malloc(Nspace * sizeof(double));
      line->Rji = (double *) malloc(Nspace * sizeof(double));

      /* --- Initialize the mutex lock for the radiative rates if there
             is more than one thread --                -------------- */

      if (input.Nthreads > 1) {
	if ((status = pthread_mutex_init(&line->rate_lock, NULL))) {
	  sprintf(messageStr, "Unable to initialize mutex_lock, status = %d",
		  status);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}
      }

      if (atmos.Stokes) {
	if (line->g_Lande_eff != 0.0 ||
	    (determinate(atom->label[i], atom->g[i], &nq, &S, &L, &Jl) &&
	     determinate(atom->label[j], atom->g[j], &nq, &S, &L, &Ju) &&
	     fabs(Ju - Jl) <= 1.0)) {

          if (line->Ncomponent > 1) {
	    sprintf(messageStr,
		    "Line %3d -> %3d: cannot treat composite line "
                    "with polarization", j, i);
	    Error(ERROR_LEVEL_2, routineName, messageStr);
	  }
	  line->polarizable = TRUE;
	} else {
	  sprintf(messageStr,
		  " -- Treating line %3d -> %3d without polarization%s\n",
		  j, i, (kr == Nline-1) ? "\n" : "");
	  Error(MESSAGE, routineName, messageStr);
	  line->polarizable = FALSE;
	}
      }
    }
  }
  /* --- Go through the bound-free transitions --      -------------- */

  atom->continuum =
    (AtomicContinuum *) malloc(Ncont * sizeof(AtomicContinuum));
  for (kr = 0;  kr < Ncont;  kr++) {
    continuum = atom->continuum + kr;
    initAtomicContinuum(continuum);
    continuum->atom = atom;
    continuum->isotope_frac = 1.0;

    getLine(atom->fp_input, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
    Nread = sscanf(inputLine, "%d %d %lf %d %s %lf",
		   &j, &i, &continuum->alpha0, &continuum->Nlambda,
		   nuDepStr, &lambdamin);
    checkNread(Nread, Nrequired=6, routineName, checkPoint=5);

    continuum->j = MAX(i, j);  continuum->i = MIN(i, j);
    j = continuum->j;
    i = continuum->i;

    lambda0 = (HPLANCK * CLIGHT)/(atom->E[j] - atom->E[i]);
    continuum->lambda0 = lambda0 / NM_TO_M;
    continuum->lambda  =
      (double *) malloc(continuum->Nlambda * sizeof(double));
    continuum->alpha   =
      (double *) malloc(continuum->Nlambda * sizeof(double));

    if (strstr(nuDepStr, "EXPLICIT")) {
      continuum->hydrogenic = FALSE;
      for (la = continuum->Nlambda-1;  la >= 0;  la--) {
	getLine(atom->fp_input, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
	Nread = sscanf(inputLine, "%lf %lf",
		        &continuum->lambda[la], &continuum->alpha[la]);
	checkNread(Nread, Nrequired=2, routineName, checkPoint=6);
      }
      for (la = 1;  la < continuum->Nlambda;  la++) {
	if (continuum->lambda[la] < continuum->lambda[la-1]) {
	  sprintf(messageStr, "Wavelength for continuum %d - %d"
		  " is not monotonous", i, j);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}
      }
    } else if (strstr(nuDepStr, "HYDROGENIC")) {
      continuum->hydrogenic = TRUE;
      if (lambdamin >= continuum->lambda0) {
	sprintf(messageStr, "Minimum wavelength for continuum %d -> %d"
		" is larger than continuum edge at: %f [nm]",
		j, i, lambdamin);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      getLambdaCont(continuum, lambdamin);
    } else {
      sprintf(messageStr,
	      "Invalid value for wavelength dependence string: %s",
	      nuDepStr);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    if (atom->active) {

      /* --- Allocate space for up- and downward radiative rates -- - */

      continuum->Rij = (double *) malloc(Nspace * sizeof(double));
      continuum->Rji = (double *) malloc(Nspace * sizeof(double));

      /* --- Initialize the mutex lock for the radiative rates if there
             is more than one thread --                -------------- */

      if (input.Nthreads > 1) {
	if ((status = pthread_mutex_init(&continuum->rate_lock, NULL))) {
	  sprintf(messageStr, "Unable to initialize mutex_lock, status = %d",
		  status);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}
      }
    }
  }

  /* --- Go through fixed transitions --               -------------- */

  if (atom->Nfixed > 0) {
    atom->ft = (FixedTransition *) malloc(Nfixed * sizeof(FixedTransition));

    for (kf = 0;  kf < Nfixed;  kf++) {
      fixed = atom->ft + kf;
      fixed->atom = atom;

      getLine(atom->fp_input, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
      Nread = sscanf(inputLine, "%d %d %lf %lf %s",
		     &j, &i, &fixed->strength,
		     &fixed->Trad, optionStr);
      checkNread(Nread, Nrequired=5, routineName, checkPoint=6);

      fixed->j = MAX(i, j);  fixed->i = MIN(i, j);
      j = fixed->j;
      i = fixed->i;

      for (kr = 0;  kr < Nline;  kr++) {
	line = atom->line + kr;
	if (line->i == i  &&  line->j == j) {
	  sprintf(messageStr,
		  "Fixed transition j = %d, i = %d duplicates active line",
		  j, i);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}
      }
      for (kr = 0;  kr < Ncont;  kr++) {
	continuum = atom->continuum + kr;
	if (continuum->i == i  &&  continuum->j == j) {
	  sprintf(messageStr, "Fixed transition j = %d,  i = %d"
		  " duplicates active continuum", j, i);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}
      }
      lambda0 = (HPLANCK * CLIGHT) / (atom->E[j] - atom->E[i]);
      fixed->lambda0 = lambda0 / NM_TO_M;

      if (atom->stage[j] == atom->stage[i])
	fixed->type = FIXED_LINE;
      else 
	fixed->type = FIXED_CONTINUUM;

      if (strstr(optionStr, "TRAD_ATMOSPHERIC"))
	fixed->option = TRAD_ATMOSPHERIC;
      else if (strstr(optionStr, "TRAD_PHOTOSPHERIC"))
	fixed->option = TRAD_PHOTOSPHERIC;
      else if (strstr(optionStr, "TRAD_CHROMOSPHERIC"))
	fixed->option = TRAD_CHROMOSPHERIC;
      else {
	sprintf(messageStr, "Invalid value for radiation temperature"
		" option string: %s", optionStr);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
    }
  }

  if (atom->active) {

      atom->popsoutFile = (char *) malloc(12 * sizeof(char));
      sprintf(atom->popsoutFile, (atom->ID[1] == ' ') ?
	      "pops.%.1s.out" : "pops.%.2s.out", atom->ID);

    if (atom->Nprd > 0) {
      if (atmos.moving && !input.PRD_angle_dep) {
	sprintf(messageStr,
		"Using angle-averaged PRD in moving atmosphere for "
                "atom %2s\n", atom->ID);
	Error(WARNING, routineName, messageStr);
      }
      if (!atmos.moving && input.PRD_angle_dep) {
	sprintf(messageStr,
		"Using angle-dependent PRD in static atmosphere for "
                "atom %2s\n", atom->ID);
	Error(WARNING, routineName, messageStr);
      }
      
      /* --- Create the array to store cross redistribution lines - - */

      if (input.XRD) {
	for (kr = 0;  kr < Nline;  kr++) {
	  line = &atom->line[kr];
	  line->xrd = (AtomicLine **) malloc(Nline * sizeof(AtomicLine *)); 
	  if (line->PRD) {
	    for (krp = 0;  krp < Nline;  krp++) {
	      if (atom->line[krp].PRD  &&
		  line->j == atom->line[krp].j  && 
		  line->i != atom->line[krp].i) {
		line->xrd[line->Nxrd] = &atom->line[krp];
		line->Nxrd++;

                /* --- Make sure wavelength quadratures match in
		       case of cross-redistribution -- -------------- */

		if (krp > kr) {
		  waveratio = atom->line[krp].lambda0 / line->lambda0;
		  atom->line[krp].qwing = waveratio * line->qwing;
		}
	      }
	    }
	    line->xrd =
	      realloc (line->xrd, line->Nxrd * sizeof(AtomicLine *));

	    sprintf(messageStr,
		    "Found %d subordinate PRD lines for line %d-%d of "
		    "atom %2s\n", line->Nxrd, line->j, line->i, atom->ID);
	    Error(MESSAGE, routineName, messageStr);
	  }
	}
      }
    }
    /* --- Get wavelength quadratures and allocate space
           for populations --                          -------------- */

    for (kr = 0;  kr < Nline;  kr++) getLambda(atom->line + kr);
    atom->n = matrix_double(Nlevel, Nspace);

    /* --- Allocate space for thread dependent quantities -- -------- */

    atom->rhth = (rhthread *) malloc(input.Nthreads * sizeof(rhthread));

    /* --- Store the offset to allow pointing back to the start of the
           collisional data in the atomic input file, and allocate
           space for rate coefficients --               ------------- */

    atom->offset_coll = ftell(atom->fp_input);
    atom->C = matrix_double(SQ(Nlevel), Nspace);

  } else {

    if (atom->popsinFile  &&
	atom->initial_solution == OLD_POPULATIONS) {
 
      atom->NLTEpops = TRUE;

      /* --- Allocate memory for Non-LTE populations -- ------------- */

      atom->n = matrix_double(atom->Nlevel, atmos.Nspace);
      readPopulations(atom);

      sprintf(messageStr,
	  " --- Read input file: %s with NLTE populations for atom %s\n",
	      atom->popsinFile, atom->ID);
      Error(MESSAGE, routineName, messageStr);
    } else {
      atom->NLTEpops = FALSE;

      /* --- Save some memory by aliasing n to nstar if passive -- -- */

      atom->n = atom->nstar;
    }
    /* --- For the inactive atom we can close the input file here. For
           the active atoms we need to read the collisional data after
           the LTE populations have been calculated (after the electron
           density has been calculated if necessary). This is done
           in routine SetLTEQuantities in ltepops.c -- -------------- */

    fclose(atom->fp_input);
  }

  sprintf(labelStr, "Read %s %2s",
	  (atom->active) ? "Active" : "Atom", atom->ID);
  getCPU(3, TIME_POLL, labelStr);
} 
/* ------- end ---------------------------- readAtom.c -------------- */

/* ------- begin -------------------------- initAtom.c -------------- */

void initAtom(Atom *atom)
{
  /* --- Initialize atomic data structure. Set pointers to NULL -- -- */

  atom->label = NULL;
  atom->popsinFile = NULL;
  atom->popsoutFile = NULL;
  atom->initial_solution = UNKNOWN;
  atom->active = FALSE;
  atom->NLTEpops = FALSE;
  atom->Nfixed = atom->Ncont = atom->Nline = atom->Nlevel = 0;
  atom->Nprd = 0;
  atom->stage = NULL;
  atom->periodic_table = 0;
  atom->abundance = atom->weight = 0.0;
  atom->g = atom->E = atom->vbroad = NULL;
  atom->C = NULL;
  atom->n = atom->nstar = NULL;
  atom->ntotal = NULL;
  atom->Gamma = NULL;
  atom->line = NULL;
  atom->continuum = NULL;
  atom->ft = NULL;
  atom->Ng_n = NULL;
  atom->fp_input = NULL;
}
/* ------- end ---------------------------- initAtom.c -------------- */

/* ------- begin -------------------------- initAtomicLine.c -------- */

void initAtomicLine(AtomicLine *line)
{
  register int n;

  /* --- Initialize transition data structure. Set pointers to NULL - */

  line->symmetric = TRUE;
  line->polarizable = FALSE;
  line->Voigt = TRUE;
  line->PRD = FALSE;
  line->vdWaals = UNSOLD;
  line->i = line->j = line->Nlambda = line->Nblue = 0;
  line->Ncomponent = 1;
  line->Nxrd = 0;
  line->lambda0 = line->isotope_frac = line->g_Lande_eff = 0.0;
  line->lambda = NULL;
  line->Aji = line->Bji = line->Bij = 0.0;
  line->Rij = line->Rji = NULL;
  line->phi = line->phi_Q = line->phi_U = line->phi_V = NULL;
  line->psi_Q = line->psi_U = line->psi_V = NULL;
  line->wphi = line->Qelast = NULL;
  line->Grad = line->cStark = 0.0;
  for (n = 0;  n < 4;  n++) line->cvdWaals[n] = 0.0;
  line->qwing = line->qcore = 0.0;
  line->c_shift = line->c_fraction = NULL;
  line->rho_prd = NULL;
  line->fp_GII = NULL;
  line->Ng_prd = NULL;
  line->atom = NULL;
  line->xrd = NULL;
}
/* ------- end ---------------------------- initAtomicLine.c -------- */

/* ------- begin -------------------------- initAtomicContinuum.c --- */

void initAtomicContinuum(AtomicContinuum *continuum)
{
  /* --- Initialize transition data structure. Set pointers to NULL - */

  continuum->hydrogenic = TRUE;
  continuum->i = continuum->j = 0;
  continuum->Nlambda = continuum->Nblue = 0;
  continuum->lambda0 = continuum->isotope_frac = 0.0;
  continuum->lambda = NULL;
  continuum->alpha0 = 0.0;
  continuum->alpha = continuum->Rij = continuum->Rji = NULL;
  continuum->atom = NULL;
}
/* ------- end ---------------------------- initAtomicContinuum.c --- */

/* ------- begin -------------------------- freeAtom.c -------------- */

void freeAtom(Atom *atom)
{
  register int kr;

  /* --- Free allocated memory for atomic data structure -- --------- */

  if (atom->label != NULL)       freeMatrix((void **) atom->label);
  if (atom->popsinFile != NULL)  free(atom->popsinFile);
  if (atom->popsoutFile != NULL) free(atom->popsoutFile);
  if (atom->stage != NULL)       free(atom->stage);
  if (atom->g != NULL)           free(atom->g);
  if (atom->E != NULL)           free(atom->E);
  if (atom->C != NULL)           freeMatrix((void **) atom->C);
  if (atom->vbroad != NULL)      free(atom->vbroad);

  /* --- Be careful here because atom->n points to atom->nstar in
         the case of LTE populations (see readAtom.c) -- ------------ */ 

  if (atom->n != atom->nstar)
    if (atom->n != NULL) freeMatrix((void **) atom->n);
  if (atom->nstar != NULL)  freeMatrix((void **) atom->nstar);
  if (atom->ntotal != NULL) free(atom->ntotal);

  if (atom->Gamma != NULL)  freeMatrix((void **) atom->Gamma);

  if (atom->line != NULL) {
    for (kr = 0;  kr < atom->Nline;  kr++)
      freeAtomicLine(atom->line + kr);
    free(atom->line);
  }
  if (atom->continuum != NULL) {
    for (kr = 0;  kr < atom->Ncont;  kr++)
      freeAtomicContinuum(atom->continuum + kr);
    free(atom->continuum);
  }
  if (atom->ft != NULL) free(atom->ft);
}
/* ------- end ---------------------------- freeAtom.c -------------- */

/* ------- begin -------------------------- freeAtomicLine.c -------- */

void freeAtomicLine(AtomicLine *line)
{
  /* --- Free allocated memory for active transition structure line - */

  if (line->lambda != NULL)  free(line->lambda);
  if (line->Rij != NULL)     free(line->Rij);
  if (line->Rji != NULL)     free(line->Rji);

  if (line->phi != NULL)     freeMatrix((void **) line->phi);
  if (line->polarizable) {
    if (line->phi_Q != NULL) freeMatrix((void **) line->phi_Q);
    if (line->phi_U != NULL) freeMatrix((void **) line->phi_U);
    if (line->phi_V != NULL) freeMatrix((void **) line->phi_V);
    if (line->psi_Q != NULL) freeMatrix((void **) line->psi_Q);
    if (line->psi_U != NULL) freeMatrix((void **) line->psi_U);
    if (line->psi_V != NULL) freeMatrix((void **) line->psi_V);
  }
  if (line->c_shift != NULL)    free(line->c_shift);
  if (line->c_fraction != NULL) free(line->c_fraction);

  if (line->wphi != NULL)    free(line->wphi);
  if (line->Qelast != NULL)  free(line->Qelast);
  if (line->rho_prd != NULL) freeMatrix((void **) line->rho_prd);
  if (line->fp_GII != NULL)  free(line->fp_GII);
}
/* ------- end ---------------------------- freeAtomicLine.c -------- */

/* ------- begin -------------------------- freeAtomicContinuum.c --- */

void freeAtomicContinuum(AtomicContinuum *continuum)
{
  /* --- Free allocated memory for AtomicContinuum structure -- ----- */

  if (continuum->lambda != NULL)  free(continuum->lambda);
  if (continuum->Rij != NULL)     free(continuum->Rij);
  if (continuum->Rji != NULL)     free(continuum->Rji);

  if (continuum->alpha != NULL)   free(continuum->alpha);
}
/* ------- end ---------------------------- freeAtomicContinuum.c --- */

/* ------- begin -------------------------- readAtomicModels.c ------ */

void readAtomicModels(void)
{
  const char routineName[] = "readAtomicModels";
  register int n, m, i;

  char    filename[MAX_LINE_SIZE],
          actionKey[MAX_KEYWORD_SIZE], popsKey[MAX_KEYWORD_SIZE],
          popsFile[MAX_LINE_SIZE], inputLine[MAX_LINE_SIZE], *atomID;
  bool_t  active, exit_on_EOF;
  int     Nread, Nrequired, checkPoint;
  FILE   *fp_atoms;
  Atom   *atom;
  Element *element;

  getCPU(2, TIME_START, NULL);

  /* --- Open input file for atomic models --          -------------- */

  if ((fp_atoms = fopen(input.atoms_input, "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s",
	    input.atoms_input);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* --- Get the number of atomic models to be read -- -------------- */

  getLine(fp_atoms, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(inputLine, "%d", &atmos.Natom);
  checkNread(Nread, Nrequired=1, routineName, checkPoint=1);

  atmos.atoms = (Atom *) malloc(atmos.Natom * sizeof(Atom));

  /* --- Read atomic data for the various model atoms -- ------------ */

  for (n = 0;  n < atmos.Natom;  n++) {

    getLine(fp_atoms, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
    Nread = sscanf(inputLine, "%s %s %s %s ",
		   filename, actionKey, popsKey, popsFile);
    checkNread(Nread, Nrequired=3, routineName, checkPoint=2);

    atomID = getAtomID(filename);
    if (n ==  0  &&  !strstr(atomID, "H ")) {
      sprintf(messageStr, "First atomic model is not hydrogen: %s",
	      atomID);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    /* --- Check for duplicate atomic ID --            -------------- */

    for (m = 0;  m < n;  m++) {
      if (strstr(atomID, atmos.atoms[m].ID)) {
	sprintf(messageStr,
		"Aready read atomic model for element %s\n", atomID);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
    }
    /* --- Set active flag. Active set to TRUE means atom will be
           treated in Non-LTE --                       -------------- */

    atom = &atmos.atoms[n];
    readAtom(atom, filename,
	     active=(strstr(actionKey, "ACTIVE") ? TRUE : FALSE));


    /* --- Set flag for initial soltion --             -------------- */

    if (strstr(popsKey, "OLD_POPULATIONS")) {
      atom->initial_solution = OLD_POPULATIONS;
    } else if (strstr(popsKey, "LTE_POPULATIONS")) {
      atom->initial_solution = LTE_POPULATIONS;
    } else if (strstr(popsKey, "ZERO_RADIATION")) {
      atom->initial_solution = ZERO_RADIATION;
    }
    /* --- If popsKey is not recognized --             -------------- */

    if (atom->initial_solution == UNKNOWN) {
      sprintf(messageStr,
	      "Unknown initial solution specified for atom: %s\n",
	      atomID);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    } 

    /* --- If input.startJ == OLD_J then enforce OLD_POPULATIONAS - - */

    if (atom->active  &&  input.startJ == OLD_J)
      atom->initial_solution = OLD_POPULATIONS;

    /* --- Copy filename for old population numbers -- -------------- */

    if (atom->initial_solution == OLD_POPULATIONS) {
      if (Nread < 4) {
	sprintf(messageStr,
		"No file with OLD_POPULATIONS specified for atom: %s\n",
		atomID);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      atom->popsinFile =
	(char *) malloc((strlen(popsFile) + 1) * sizeof(char));
      strcpy(atom->popsinFile, popsFile);
    }
  }
  fclose(fp_atoms);  

  /* --- Create an alias to the hydrogen atom structure -- ---------- */

  atmos.H = &atmos.atoms[0];
  
  /* --- Store pointers to models in element structures -- ---------- */

  for (n = 0;  n < atmos.Natom;  n++)
    atmos.elements[atmos.atoms[n].periodic_table].model = &atmos.atoms[n];

  /* --- Redistribute the hydrogen populations read in with the
         atmosphere over the atmospheric Hydrogen model -- ---------- */

  distribute_nH();
  getCPU(2, TIME_POLL, "Read atomic input");
}

/* ------- end ---------------------------- readAtomicModels.c ------ */

/* ------- begin -------------------------- getAtomID.c ------------- */

char *getAtomID(char *atom_file)
{
  const char routineName[] = "getAtomID";
  register int n;
  static char atomID[ATOM_ID_WIDTH + 1];

  char   inputLine[MAX_LINE_SIZE];
  bool_t exit_on_EOF;
  FILE  *fp_atom;

  if ((fp_atom = fopen(atom_file, "r")) == NULL) {
    sprintf(messageStr, "Unable to open inputfile %s", atom_file);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* --- Read atom ID and convert to uppercase --     -------------- */
 
  getLine(fp_atom, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  sscanf(inputLine, "%2s", atomID);

  for (n = 0;  n < (int) strlen(atomID);  n++)
    atomID[n] = toupper(atomID[n]);
  if (strlen(atomID) == 1) strcat(atomID, " ");

  fclose(fp_atom);
  return atomID;
}
/* ------- end ---------------------------- getAtomID.c ------------- */
