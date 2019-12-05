/* ------- file: -------------------------- rawAtom.c ---------------

       Version:       rh1.0, tools
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Jun 19 10:13:25 2001 --

       --------------------------                      ----------RH-- */

/* --- Reads atomic model into Atom structure pointed to by atom. - - */

 
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "atomweights.h"

#define  COMMENT_CHAR        "#"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

char messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- rawAtom.c --------------- */

void rawAtom(Atom *atom, char *atomFileName)
{
  const char routineName[] = "rawAtom";
  register int kr, kf, la, n;

  char    inputLine[MAX_LINE_SIZE], shapeStr[20], vdWstr[20], nuDepStr[20],
          errorStr[80], symmStr[20], optionStr[20], *match;
  bool_t  exit_on_EOF;
  int     i, j, Nlevel, Nrad, Nline, Ncont, Nfixed, Nread, Nrequired,
          checkPoint;
  double  f, C, lambda0, lambdamin, dlamb;
  FILE   *atomFile;
  AtomicLine *line;
  AtomicContinuum *continuum;
  FixedTransition *fixed;

  /* --- Open the data file for current model atom --  -------------- */

  if ((atomFile = fopen(atomFileName, "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", atomFileName);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* --- Read atom ID and convert to uppercase --      -------------- */
 
  getLine(atomFile, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(inputLine, "%2s", atom->ID);
  checkNread(Nread, Nrequired=1, routineName, checkPoint=1);
  for (n = 0;  n < (int) strlen(atom->ID);  n++)
    atom->ID[n] = toupper(atom->ID[n]);
  if (strlen(atom->ID) == 1) strcat(atom->ID, " ");
 
  /* --- Get the atomic weight --                      -------------- */

  for (n = 0;  n < sizeof(atomweight)/sizeof(struct AtomWeight);  n++) {
    if ((match = strstr(atom->ID, atomweight[n].ID))) {
      atom->weight = atomweight[n].weight;
      break;
    }
  }
  if (!match) {
    sprintf(messageStr, "Found no matching element for %s", atom->ID);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* --- Get Number of levels, lines fixed transitions, and continua  */
 
  getLine(atomFile, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(inputLine, "%d %d %d %d",
		 &atom->Nlevel, &atom->Nline, &atom->Ncont, &atom->Nfixed);
  checkNread(Nread, Nrequired=4, routineName, checkPoint=2);
  Nlevel = atom->Nlevel;
  Nline  = atom->Nline;  Ncont  = atom->Ncont;  Nfixed = atom->Nfixed;
  Nrad   = Nline + Ncont;

  atom->E = (double *) malloc(Nlevel * sizeof(double));
  atom->g = (double *) malloc(Nlevel * sizeof(double));
  atom->stage = (int *)   malloc(Nlevel * sizeof(int));
  atom->label = (char **) malloc(Nlevel * sizeof(char *));

  /* --- Read in the level energies, statistical weights, labels,
         and ionization stage --                       -------------- */

  for (i = 0;  i < Nlevel;  i++) {
    atom->label[i] = (char *) calloc((ATOM_LABEL_WIDTH+1), sizeof(char)); 
    getLine(atomFile, COMMENT_CHAR, inputLine , exit_on_EOF=TRUE);
    Nread = sscanf(inputLine, "%lf %lf '%20c' %d",
      &atom->E[i], &atom->g[i], atom->label[i], &atom->stage[i]);
    checkNread(Nread, Nrequired=4, routineName, checkPoint=3);

    atom->E[i] *= (HPLANCK * CLIGHT) / CM_TO_M;
  }
  if (atom->stage[Nlevel-1] != (atom->stage[Nlevel-2] + 1)) {
    sprintf(messageStr, "Found no overlying continuum for atom %s", atom->ID);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  C = 2 * PI * (Q_ELECTRON/EPSILON_0) * (Q_ELECTRON/M_ELECTRON) / CLIGHT;

  /* --- Go through the bound-bound transitions --     -------------- */

  atom->Nprd = 0;
  atom->line = (AtomicLine *) malloc(Nline * sizeof(AtomicLine));
  for (kr = 0;  kr < Nline;  kr++) {
    line = atom->line + kr;
    line->atom = atom;
    line->isotope_frac = 1.0;

    getLine(atomFile, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
    Nread = sscanf(inputLine,
		   "%d %d %lf %s %d %s %lf %lf %s %lf %lf %lf %lf %lf %lf %lf",
		   &j, &i, &f, shapeStr,
		   &line->Nlambda, symmStr, &line->qcore,
		   &line->qwing, vdWstr,
		   &line->cvdWaals[0], &line->cvdWaals[1],
		   &line->cvdWaals[2], &line->cvdWaals[3],
		   &line->Grad, &line->cStark, &line->g_Lande_eff);
    checkNread(Nread, Nrequired=15, routineName, checkPoint=4);
    if (Nread == 15) line->g_Lande_eff = 0.0;

    line->j = MAX(i, j);  line->i = MIN(i, j);
    i = line->i;
    j = line->j;

    if (!strstr(shapeStr, "PRD") &&
	!strstr(shapeStr, "VOIGT") && !strstr(shapeStr, "GAUSS")) {
      sprintf(messageStr, "Invalid value for line-shape string: %s",
	      shapeStr);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    if (strstr(shapeStr, "PRD")) {
      atom->Nprd++;
      line->PRD = TRUE;
    }
    if (strstr(shapeStr, "GAUSS")) line->Voigt = FALSE;

    lambda0 = (HPLANCK * CLIGHT) / (atom->E[j] - atom->E[i]);
    line->Aji = C / SQ(lambda0) * (atom->g[i] / atom->g[j]) * f;
    line->Bji = CUBE(lambda0) / (2.0 * HPLANCK * CLIGHT) * line->Aji;
    line->Bij = (atom->g[j] / atom->g[i]) * line->Bji;
    line->lambda0 = lambda0 / NM_TO_M;

    if (strstr(vdWstr, "PARAMTR"))
      line->vdWaals = RIDDER_RENSBERGEN;
    else {
      line->vdWaals = UNSOLD;
      line->cvdWaals[3] = line->cvdWaals[1] = 0.0;
    }
    line->symmetric = (strstr(symmStr, "ASYMM")) ? FALSE : TRUE;
  }
  /* --- Go through the bound-free transitions --      -------------- */


  atom->continuum =
    (AtomicContinuum *) malloc(Ncont * sizeof(AtomicContinuum));
  for (kr = 0;  kr < Ncont;  kr++) {
    continuum = atom->continuum + kr;
    continuum->atom = atom;
    continuum->isotope_frac = 1.0;

    getLine(atomFile, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
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
	getLine(atomFile, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
	Nread = sscanf(inputLine, "%lf %lf",
		        &continuum->lambda[la], &continuum->alpha[la]);
	checkNread(Nread, Nrequired=2, routineName, checkPoint=6);
      }
      for (la = 1;  la < continuum->Nlambda;  la++) {
	if (continuum->lambda[la] < continuum->lambda[la-1]) {
	  sprintf(messageStr,
		  "Continuum %d does not have monotonous wavelengths",
		  kr - Nline);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}
      }
    } else {
      continuum->hydrogenic = TRUE;
      if (lambdamin >= continuum->lambda0) {
	sprintf(messageStr,
		"Minimum wavelength for continuum %d too long", kr - Nline);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      dlamb   = (continuum->lambda0 - lambdamin) / (continuum->Nlambda - 1);
      continuum->lambda[0] = lambdamin;
      for (la = 1;  la < continuum->Nlambda;  la++)
	continuum->lambda[la] = continuum->lambda[la-1] + dlamb;
    }
  }
  /* --- Go through fixed transitions --               -------------- */

  if (atom->Nfixed > 0) {
    atom->ft = (FixedTransition *) malloc(Nfixed * sizeof(FixedTransition));

    for (kf = 0;  kf < Nfixed;  kf++) {
      fixed = atom->ft + kf;

      getLine(atomFile, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
      Nread = sscanf(inputLine, "%d %d %lf %lf %s",
		     &j, &i, &fixed->strength,
		     &fixed->Trad, optionStr);
      checkNread(Nread, Nrequired=5, routineName, checkPoint=7);

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
	sprintf(messageStr, "Inavlid value for TRAD option: %s", optionStr);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
    }
  }

  fclose(atomFile);
} 
/* ------- end ---------------------------- rawAtom.c --------------- */
