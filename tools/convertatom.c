/* ------- file: -------------------------- convertatom.c -----------

       Version:       rh1.0, tools
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Apr 29 03:46:50 2011 --

       --------------------------                      ----------RH-- */


/* --- Convert atom file in MULTI format to rh format.

       To preview:     convertatom inputFile
       To store:       convertatom inputFile > outputFile,
              or       convertatom inputFile outputFile

       --                                              -------------- */ 
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"

#define MULTI_COMMENT_CHAR   "*"


/* --- Function prototypes --                          -------------- */

void MULTIatom(FILE *atomFile, struct Atom *atom);
void writeModelAtom(struct Atom *atom, FILE *fpOut);
void writeCollisions(FILE *atomFile, struct Atom *atom, FILE *fpOut);


/* --- Global variables --                             -------------- */

CommandLine commandline;
char messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- convertAtom.c ----------- */

int main(int argc, char *argv[])
{
  Atom atom;
  FILE   *fpOut, *atomFile;

  commandline.quiet = FALSE;
  commandline.logfile = stderr;

  if (argc >= 3)
    fpOut = fopen(argv[2], "w");
  else if (argc == 2)
    fpOut = stdout;
  else {
    fprintf(stderr, "Usage: %s inFile [outFile]\n", argv[0]);
    exit(0);
  }
  /* --- Open the data file for model atom --          -------------- */

  if ((atomFile = fopen(argv[1], "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", argv[1]);
    Error(ERROR_LEVEL_2, "convertAtom", messageStr);
  }

  MULTIatom(atomFile, &atom);
  writeModelAtom(&atom, fpOut);
  writeCollisions(atomFile, &atom, fpOut);

  fclose(atomFile);  fclose(fpOut);
}
/* ------- end ---------------------------- convertAtom.c ----------- */

/* ------- begin -------------------------- MULTIatom.c ------------- */

/* Reads atomic model in MULTI format and stores it in atomic model
   structure atom.

 * There are many differences between the two formats, but fortunately
   there are enough quantities similar to make automatic conversion
   feasible.

 * - Units in the new code are in SI, while MULTI's input is mostly
     in cgs units.
   - Index counts start with 0 in the new code (as is customary
     in IDL and C) and not with 1 (as in FORTRAN).

 * You probably will have to check the output and change some things
   manually.

 * Han Uitenbroek
 * Last modified: Nov 5, 1995
 */
 
void MULTIatom(FILE *atomFile, Atom *atom)
{
  register int n, kr, kf, la, krp;

  char   input[MAX_LINE_SIZE];
  bool_t exit_on_EOF;
  int    i, j, Nrad, Nread, phot, Trad_option, iw, lap, Nbb;
  double C, f, A0, Trad, lambda0, lambdaMin;
  AtomicLine *line;
  AtomicContinuum *continuum;
  FixedTransition *fixed;

  /* --- Read atom ID --                               -------------- */

  getLine(atomFile, MULTI_COMMENT_CHAR, input, exit_on_EOF=TRUE);
  Nread = sscanf(input, "%2s", atom->ID);
  if (strlen(atom->ID) == 1) strcat(atom->ID, " ");
  for (n = 0;  n < ATOM_ID_WIDTH;  n++)
    atom->ID[n] = toupper(atom->ID[n]);

  /* --- Read atomic abundance and weight --           -------------- */

  getLine(atomFile, MULTI_COMMENT_CHAR, input, exit_on_EOF=TRUE);
  Nread = sscanf(input, "%lf %lf", &atom->abundance, &atom->weight);
  atom->abundance = POW10(atom->abundance - 12.0);

  /* --- Get Number of levels, lines fixed transitions, and continua  */

  getLine(atomFile, MULTI_COMMENT_CHAR, input, exit_on_EOF=TRUE);
  Nread = sscanf(input, "%d %d %d %d",
		 &atom->Nlevel, &atom->Nline, &atom->Ncont, &atom->Nfixed);

  Nrad = atom->Nline + atom->Ncont;

  /* --- Read in the level energies, statistical weights, labels,
         and ionization stage --                       -------------- */

  atom->E = (double *) malloc(atom->Nlevel * sizeof(double));
  atom->g = (double *) malloc(atom->Nlevel * sizeof(double));
  atom->label = (char **) malloc(atom->Nlevel * sizeof(char *));
  atom->stage = (int *) malloc(atom->Nlevel * sizeof(int));

  for (i = 0;  i < atom->Nlevel;  i++) {
    atom->label[i] = (char *) malloc((ATOM_LABEL_WIDTH+1) * sizeof(char)); 
    getLine(atomFile, MULTI_COMMENT_CHAR, input , exit_on_EOF=TRUE);
    Nread = sscanf(input, "%lf %lf '%20c' %d",
      &atom->E[i], &atom->g[i], atom->label[i], &atom->stage[i]);

    *(atom->label[i] + ATOM_LABEL_WIDTH) = '\0';
    atom->E[i] *= (HPLANCK * CLIGHT) / CM_TO_M;
    atom->stage[i]--;
  }
  C = 2 * PI * (Q_ELECTRON/EPSILON_0) * (Q_ELECTRON/M_ELECTRON) / CLIGHT;

  /* --- Go through the bound-bound transitions --     -------------- */

  atom->line = (AtomicLine *) malloc(atom->Nline * sizeof(AtomicLine));
  for (kr = 0;  kr < atom->Nline;  kr++) {
    line = atom->line + kr;

    line->Voigt   = TRUE;
    line->PRD     = FALSE;
    line->vdWaals = UNSOLD;

    getLine(atomFile, MULTI_COMMENT_CHAR, input, exit_on_EOF=TRUE);
    Nread = sscanf(input, "%d %d %lf %d %lf %lf %d %lf %lf %lf",
		   &j, &i, &f, &line->Nlambda, &line->qwing,
		   &line->qcore, &iw,
		   &line->Grad, &line->cvdWaals[0], &line->cStark);
    
    line->j = MAX(i, j);  line->i = MIN(i, j);
    j = --line->j;  i = --line->i;

    lambda0 = (HPLANCK * CLIGHT) / (atom->E[j] - atom->E[i]);
    line->Aji = C / SQ(lambda0) * (atom->g[i] / atom->g[j]) * f;
    line->lambda0 = lambda0 / NM_TO_M;

    if ((line->qcore < 0.0) || (line->qwing < 0.0)) {
      line->lambda = (double *) malloc(line->Nlambda * sizeof(double));
      for (la = 0;  la < line->Nlambda;  la++)
	Nread = fscanf(atomFile, "%lf", line->lambda + la);
      line->symmetric = FALSE;
      line->qcore = fabs(line->qcore);
      line->qwing = fabs(line->qwing);
    } else {
      line->lambda = NULL;
      line->symmetric = TRUE;
    }
    line->cvdWaals[2] = line->cvdWaals[0];
    line->cStark *= -CUBE(CM_TO_M);
  }
  /* --- Go through the bound-free transitions --      -------------- */

  atom->continuum =
    (AtomicContinuum *) malloc(atom->Ncont * sizeof(AtomicContinuum));
  for (kr = 0;  kr < atom->Ncont;  kr++) {
    continuum = atom->continuum + kr;

    getLine(atomFile, MULTI_COMMENT_CHAR, input, exit_on_EOF=TRUE);
    Nread = sscanf(input, "%d %d %lf %d %lf",
		   &j, &i, &continuum->alpha0,
		   &continuum->Nlambda, &lambdaMin);

    continuum->j = MAX(i, j);  continuum->i = MIN(i, j);
    j = --continuum->j;  i = --continuum->i;

    lambda0 = (HPLANCK * CLIGHT)/(atom->E[j] - atom->E[i]);
    continuum->lambda0 = lambda0 / NM_TO_M;
    continuum->alpha0 *= SQ(CM_TO_M);

    if (lambdaMin < 0.0) {
      continuum->hydrogenic = FALSE;
      continuum->lambda =
	(double *) malloc(continuum->Nlambda * sizeof(double));
      continuum->alpha  =
	(double *) malloc(continuum->Nlambda * sizeof(double));
      for (la = 0;  la < continuum->Nlambda;  la++) {
	lap = continuum->Nlambda - la - 1;
	getLine(atomFile, MULTI_COMMENT_CHAR, input, exit_on_EOF=TRUE);
	sscanf(input, "%lf %lf", continuum->lambda+lap, continuum->alpha+lap);
	continuum->lambda[lap] *= 0.1;
        continuum->alpha[lap]  *= SQ(CM_TO_M);
      } 
    } else {
      continuum->hydrogenic = TRUE;
      continuum->alpha = NULL;
      continuum->lambda = (double *) malloc(sizeof(double));
      continuum->lambda[0] = lambdaMin / 10.0;
    }
  }
  /* --- Go through the fixed transitions, which may be either
         bound-bound (phot == 0) or bound-free (phot == 1) -- ------- */

  atom->ft =
    (FixedTransition *) malloc(atom->Nfixed * sizeof(FixedTransition));
  for (kf = 0;  kf < atom->Nfixed;  kf++) {
    fixed = atom->ft + kf;

    getLine(atomFile, MULTI_COMMENT_CHAR, input, exit_on_EOF=TRUE);
    Nread = sscanf(input, "%d %d %d %lf %lf %d",
		   &j, &i, &phot, &A0, &fixed->Trad, &Trad_option);

    fixed->j = MAX(i, j);  fixed->i = MIN(i, j);
    j = --fixed->j;  i = --fixed->i;
    
    lambda0 = (HPLANCK * CLIGHT)/(atom->E[j] - atom->E[i]);
    fixed->lambda0 = lambda0 / NM_TO_M;
    if (phot) {
      fixed->type = FIXED_LINE;
      fixed->strength = A0 * SQ(CM_TO_M);
    } else {
      fixed->type = FIXED_CONTINUUM;
      fixed->strength = A0;
    }
    switch (Trad_option) {
    case 1: fixed->option = TRAD_ATMOSPHERIC;    break;
    case 2: fixed->option = TRAD_PHOTOSPHERIC;   break;
    case 3: fixed->option = TRAD_CHROMOSPHERIC;  break;
    }
  }
} 
/* ------- end ---------------------------- MULTIatom.c ------------- */
