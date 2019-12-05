/* ------- file: -------------------------- modelatom.c -------------

       Version:       rh1.0, tools
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Feb  3 09:55:09 2000 --

       --------------------------                      ----------RH-- */

/* --- Writes model atom from Atom structure atom. --  -------------- */

 
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"

#define COMMENT_CHAR         "#"
#define MULTI_COMMENT_CHAR   "*"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

char messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- writeModelAtom.c -------- */

void writeModelAtom(Atom *atom, FILE *fp_out)
{
  register int kr, kf, la, krp;

  char   shapeStr[20], symmStr[6], vdwStr[] = "UNSOLD", optionStr[20];
  int    i, j;
  double abundance, Ecm, f, C;
  AtomicLine *line, *another_line;
  AtomicContinuum *continuum;
  FixedTransition *fixed;

  fprintf(fp_out, "%1s <<One-line atomic model description>>\n", COMMENT_CHAR);
  fprintf(fp_out, "  %2s\n\n", atom->ID);

  abundance = log10(atom->abundance) + 12;
  fprintf(stderr,
	  "\n Warning: abundance value %7.4f and weight %7.3f are ignored\n\n",
	  abundance, atom->weight);
  fprintf(fp_out, "%1s Nlevel  Nline   Ncont   Nfixed\n", COMMENT_CHAR);
  fprintf(fp_out, "  %3d    %4d     %3d       %2d\n\n",
	  atom->Nlevel, atom->Nline, atom->Ncont, atom->Nfixed);

  fprintf(fp_out,
	  "%1s   E[cm^-1]    g           label[20]         stage   levelNo\n",
          COMMENT_CHAR);
  fprintf(fp_out, "%1s                     '|----|----|----|----'\n",
          COMMENT_CHAR);
  for (i = 0;  i < atom->Nlevel;  i++) {
    Ecm = CM_TO_M * atom->E[i] / (HPLANCK * CLIGHT);
    fprintf(fp_out, " %10.3f  %6.2f   '%-20s'    %1d     %3d\n",
	    Ecm, atom->g[i], atom->label[i], atom->stage[i], i);
  }

  C = 2 * PI * (Q_ELECTRON/EPSILON_0) * (Q_ELECTRON/M_ELECTRON) / CLIGHT;

  fprintf(fp_out,
	  "\n%1s  j    i       f     type  Nlambda symmetr  qcore  qwing",
          COMMENT_CHAR);
  fprintf(fp_out, "  vdWapprx          vdWaals            radiative    Stark");
  fprintf(fp_out,
	  "\n%1s                                                      ",
          COMMENT_CHAR);
  fprintf(fp_out, "                  H            He\n");
  for (kr = 0;  kr < atom->Nline;  kr++) {
    line = atom->line + kr;
    i = line->i;  j = line->j;
    if (!line->PRD) {
      if (line->Voigt)
	strcpy(shapeStr, "VOIGT");
      else
	strcpy(shapeStr, "GAUSS");
    } else
      strcpy(shapeStr, " PRD ");

    if (line->symmetric)
      strcpy(symmStr, " SYMM");
    else
      strcpy(symmStr, "ASYMM");

    if (line->Grad == 0.0) {
      for (krp = 0;  krp < atom->Nline;  krp++) {
	if ((another_line = atom->line + krp)->j == j)
	  line->Grad += another_line->Aji;
      }
      for (kf = 0;  kf < atom->Nfixed;  kf++) {
        fixed = atom->ft + kf;
	if ((fixed->type == ATOMIC_LINE) && (fixed->j == j)) {
	  line->Grad += C / SQ(fixed->lambda0) *
	    (atom->g[fixed->i] / atom->g[fixed->j]) * fixed->strength;
	}
      }
    }
    f = line->Aji * (atom->g[j] / atom->g[i]) *
      SQ(line->lambda0 * NM_TO_M) / C;
    fprintf(fp_out, " %3d  %3d  %9.3E  %s    %3d    %5s  %4.1f  %6.1f",
	    j, i, f, shapeStr, line->Nlambda, symmStr,
	    line->qcore, line->qwing);
    fprintf(fp_out, "   %-7s  %5.3f  %5.3f  %5.3f  %5.3f  %8.2E  %9.2E\n",
	    vdwStr, line->cvdWaals[0], line->cvdWaals[1],
	    line->cvdWaals[2], line->cvdWaals[3], line->Grad, line->cStark);
  }

  fprintf(fp_out, "\n%1s Photoionization Cross Sectionss\n\n", COMMENT_CHAR);
  fprintf(fp_out,
     "%1s j   i  alpha [m^-2]     Nlambda     Wavel. Dep.   lamb_min [nm]",
     COMMENT_CHAR);
  fprintf(fp_out, "\n%1s\n", COMMENT_CHAR);
  for (kr = 0;  kr < atom->Ncont;  kr++) {
    continuum = atom->continuum + kr;
    i = continuum->i;
    j = continuum->j;
    if (continuum->hydrogenic)
      strcpy(shapeStr, "HYDROGENIC");
    else
      strcpy(shapeStr, " EXPLICIT ");

    fprintf(fp_out, "\n%1s   %s\n", COMMENT_CHAR, atom->label[i]);
    fprintf(fp_out, " %2d  %2d     %10.3E      %3d       %s     %8.3f\n",
	    j, i, continuum->alpha0, continuum->Nlambda, shapeStr,
	    continuum->lambda[0]);

    if (!continuum->hydrogenic) {
      for (la = continuum->Nlambda-1;  la >= 0;  la--) {
	fprintf(fp_out, " %8.3f   %10.3E\n",
		continuum->lambda[la], continuum->alpha[la]);
      }
    }
  }

  fprintf(fp_out, "\n%1s Fixed Transitions\n", COMMENT_CHAR);
  fprintf(fp_out, "\n%1s j   i      Strength       Trad       Option\n\n",
	  COMMENT_CHAR);

  for (kf = 0;  kf < atom->Nfixed;  kf++) {
    fixed = atom->ft + kf;
    switch (fixed->option) {
    case TRAD_ATMOSPHERIC:   strcpy(optionStr, "TRAD_ATMOSPHERIC");   break;
    case TRAD_PHOTOSPHERIC:  strcpy(optionStr, "TRAD_PHOTOSPHERIC");  break;
    case TRAD_CHROMOSPHERIC: strcpy(optionStr, "TRAD_CHROMOSPHERIC"); break;
    }
    fprintf(fp_out, " %2d  %2d     %10.3E      %6.1f     %s\n",
	    fixed->j, fixed->i, fixed->strength, fixed->Trad, optionStr);
  }
}
/* ------- end ---------------------------- writeModelAtom.c -------- */

/* ------- begin -------------------------- writeCollisions.c ------- */

void writeCollisions(FILE *atomFile, Atom *atom, FILE *fp_out)
{
  register int n;

  char    collisionRoutine[20], keyword[20], input[MAX_LINE_SIZE];
  bool_t  exit_on_EOF;
  int     nread, Ntemp, i1, i2, i, j;
  double *temp = NULL;

  /* --- Transforms collisional input for GENCOL routine in MULTI
         to new format. Currently only OHM, CI, CE, and CP entries are
         supported --                                 --------------- */

  getLine(atomFile, MULTI_COMMENT_CHAR, input, exit_on_EOF=TRUE);
  nread = sscanf(input, "%s", collisionRoutine);
  if ((strcmp(collisionRoutine, "GENCOL") != 0)) {
    fprintf(stderr, "\n Error: can only treat collisional data from\n ");
    fprintf(stderr, " GENCOL, not from routine %s \n", collisionRoutine);
    exit(0);
  }
  fprintf(fp_out, "\n# Collisional data\n");

  while((getLine(atomFile, MULTI_COMMENT_CHAR, input,
		 exit_on_EOF=FALSE) != EOF)) {
    nread = sscanf(input, "%s", keyword);

    if (strcmp(keyword, "TEMP") == 0) {
      getLine(atomFile, MULTI_COMMENT_CHAR, input, exit_on_EOF=TRUE);
      sscanf(strtok(input, " "), "%d", &Ntemp);
      temp = (double *) realloc(temp, Ntemp * sizeof(double));

      fprintf(fp_out, "\n TEMP   %2d       ", Ntemp);
      for (n = 0;  n < Ntemp;  n++) {
	sscanf(strtok(NULL, " "), "%lf", temp+n);
	fprintf(fp_out, "%9.1f   ", temp[n]);
      }
      fprintf(fp_out, "\n\n");

    } else if (strcmp(keyword, "OHM") == 0) { 
      getLine(atomFile, MULTI_COMMENT_CHAR, input, exit_on_EOF=TRUE);
      sscanf(strtok(input, " "), "%d", &i1);
      sscanf(strtok(NULL, " "), "%d", &i2);

      j = MAX(i1, i2) - 1;
      i = MIN(i1, i2) - 1;
      fprintf(fp_out, " OMEGA  %2d  %2d   ", j, i);
      for (n = 0;  n < Ntemp;  n++) {
	sscanf(strtok(NULL, " "), "%lf", temp+n);
	fprintf(fp_out, "%9.3E   ", temp[n]);
      }
      fprintf(fp_out, "\n");

    } else if (strcmp(keyword, "CE") == 0) {
      getLine(atomFile, MULTI_COMMENT_CHAR, input, exit_on_EOF=TRUE);
      sscanf(strtok(input, " "), "%d", &i1);
      sscanf(strtok(NULL, " "), "%d", &i2);

      j = MAX(i1, i2) - 1;
      i = MIN(i1, i2) - 1;
      fprintf(fp_out, " CE     %2d  %2d   ", j, i);
      for (n = 0;  n < Ntemp;  n++) {
	sscanf(strtok(NULL, " "), "%lf", temp+n);
	fprintf(fp_out, "%9.3E   ", temp[n] * CUBE(CM_TO_M));
      }
      fprintf(fp_out, "\n");

    } else if (strcmp(keyword, "CI") == 0) {
      getLine(atomFile, MULTI_COMMENT_CHAR, input, exit_on_EOF=TRUE);
      sscanf(strtok(input, " "), "%d", &i1);
      sscanf(strtok(NULL, " "), "%d", &i2);

      j = MAX(i1, i2) - 1;
      i = MIN(i1, i2) - 1;
      fprintf(fp_out, " CI     %2d  %2d   ", i, j);
      for (n = 0;  n < Ntemp;  n++) {
	sscanf(strtok(NULL, " "), "%lf", temp+n);
	fprintf(fp_out, "%9.3E   ", temp[n] * CUBE(CM_TO_M));
      }
      fprintf(fp_out, "\n");

    } else if (strcmp(keyword, "CP") == 0) {
      getLine(atomFile, MULTI_COMMENT_CHAR, input, exit_on_EOF=TRUE);
      sscanf(strtok(input, " "), "%d", &i1);
      sscanf(strtok(NULL, " "), "%d", &i2);

      j = MAX(i1, i2) - 1;
      i = MIN(i1, i2) - 1;
      fprintf(fp_out, " CP     %2d  %2d   ", j, i);
      for (n = 0;  n < Ntemp;  n++) {
	sscanf(strtok(NULL, " "), "%lf", temp+n);
	fprintf(fp_out, "%9.3E   ", temp[n] * CUBE(CM_TO_M));
      }
      fprintf(fp_out, "\n");

    } else if (strcmp(keyword, "CH0") == 0) {
      getLine(atomFile, MULTI_COMMENT_CHAR, input, exit_on_EOF=TRUE);
      sscanf(strtok(input, " "), "%d", &i1);
      sscanf(strtok(NULL, " "), "%d", &i2);

      j = MAX(i1, i2) - 1;
      i = MIN(i1, i2) - 1;
      fprintf(fp_out, " CH0    %2d  %2d   ", j, i);
      for (n = 0;  n < Ntemp;  n++) {
	sscanf(strtok(NULL, " "), "%lf", temp+n);
	fprintf(fp_out, "%9.3E   ", temp[n] * CUBE(CM_TO_M));
      }
      fprintf(fp_out, "\n");


    } else if (strcmp(keyword, "CH+") == 0) {
      getLine(atomFile, MULTI_COMMENT_CHAR, input, exit_on_EOF=TRUE);
      sscanf(strtok(input, " "), "%d", &i1);
      sscanf(strtok(NULL, " "), "%d", &i2);

      j = MAX(i1, i2) - 1;
      i = MIN(i1, i2) - 1;
      fprintf(fp_out, " CH+    %2d  %2d   ", i, j);
      for (n = 0;  n < Ntemp;  n++) {
	sscanf(strtok(NULL, " "), "%lf", temp+n);
	fprintf(fp_out, "%9.3E   ", temp[n] * CUBE(CM_TO_M));
      }
      fprintf(fp_out, "\n");

    } else if (strcmp(keyword, "END") == 0) {
      fprintf(fp_out, "\n END");
      break;
    } else {
      sprintf(messageStr, "Unknown keyword: %s", keyword);
      Error(ERROR_LEVEL_2, "writeCollisions", messageStr);
    }
  }
}
/* ------- end ---------------------------- writeCollisions.c ------- */
